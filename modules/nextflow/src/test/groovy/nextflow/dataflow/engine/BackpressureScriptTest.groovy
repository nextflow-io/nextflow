/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nextflow.dataflow.engine

import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.TimeUnit

import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Global
import nextflow.NF
import nextflow.NextflowMeta
import nextflow.Session
import nextflow.SysEnv
import nextflow.exception.ProcessFailedException
import nextflow.executor.Executor
import nextflow.executor.ExecutorFactory
import nextflow.executor.local.NativeTaskHandler
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.ScriptBinding
import nextflow.script.ScriptLoaderFactory
import nextflow.script.ScriptMeta
import nextflow.script.ScriptType
import nextflow.trace.WorkflowStats
import nextflow.util.Duration
import spock.lang.Specification
import spock.lang.TempDir
import spock.lang.Timeout
import test.ScriptHelper

/**
 * Runs the deadlock cases of the backpressure design with backpressure enabled, on a real
 * task polling monitor with two task slots, and checks that each run terminates and keeps
 * the process backlog bounded.
 *
 * Every run, including a failed one, must end with every process backlog back to zero and
 * every gate open, so that an accounting leak, which would hang a real run, fails the test.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Timeout(value = 60, unit = TimeUnit.SECONDS)
class BackpressureScriptTest extends Specification {

    static final int HW = 8

    static final int N = 1000

    /**
     * A task is counted both as preparing and pending between its scheduling and the end of
     * its preparation, so the sampled count may exceed the mark by the tasks being prepared.
     */
    static final int MAX_BUSY = 2 * HW

    @TempDir
    Path workDir

    def setup() {
        reset()
        SysEnv.push(NXF_DATAFLOW_HIGH_WATER: String.valueOf(HW))
    }

    private static void reset() {
        TaskProcessor.reset()
        ScriptMeta.reset()
        Global.reset()
        NF.init()
    }

    def cleanup() {
        SysEnv.pop()
    }

    static class TestExecutor extends Executor {
        @Override
        protected TaskMonitor createTaskMonitor() {
            return TaskPollingMonitor.create(getSession(), getConfig(), getName(), 2, Duration.of('5ms'))
        }

        @Override
        TaskHandler createTaskHandler(TaskRun task) {
            // an anonymous subclass, since the constructor is protected
            return new NativeTaskHandler(task, this) {}
        }
    }

    static class TestExecutorFactory extends ExecutorFactory {
        @Override
        protected Class<? extends Executor> getExecutorClass(String executorName) {
            return TestExecutor
        }

        @Override
        protected boolean isTypeSupported(ScriptType type, Object executor) {
            return true
        }
    }

    static class TestSession extends Session {
        TestSession(Map<String,Object> config) {
            super(config)
        }

        @Override
        Session start() {
            setExecutorFactory(new TestExecutorFactory())
            return super.start()
        }
    }

    /**
     * Samples the largest backlog of each process during a run, and keeps each backlog
     * since a process is deregistered when it terminates.
     */
    static class Sampler extends Thread {
        final Backpressure bp
        final Map<String,ProcessBacklog> backlogs = new ConcurrentHashMap<>()
        final Map<String,Integer> busy = new ConcurrentHashMap<>()
        final Map<String,Integer> ready = new ConcurrentHashMap<>()
        volatile boolean done

        Sampler(Backpressure bp) {
            this.bp = bp
            setDaemon(true)
        }

        @Override
        void run() {
            while( !done ) {
                for( ProcessBacklog it : bp.processes() ) {
                    backlogs.putIfAbsent(it.getName(), it)
                    busy.merge(it.getName(), it.getPreparing() + it.getPending(), { a, b -> Math.max(a, b) })
                    ready.merge(it.getName(), it.ready(100_000), { a, b -> Math.max(a, b) })
                }
                try {
                    Thread.sleep(1)
                }
                catch( InterruptedException e ) {
                    return
                }
            }
        }
    }

    static class Run {
        final Object result
        final Session session
        final Sampler sampler

        Run(Object result, Session session, Sampler sampler) {
            this.result = result
            this.session = session
            this.sampler = sampler
        }

        Backpressure backpressure() {
            return session.getBackpressure()
        }

        Object value() {
            return ((DataflowReadChannel) result).getVal()
        }

        int maxBusy(String name) {
            return sampler.busy.getOrDefault(name, 0)
        }

        int maxReady(String name) {
            return sampler.ready.getOrDefault(name, 0)
        }

        Backpressure.Gate gate(String source) {
            for( Backpressure.Gate gate : backpressure().gates() ) {
                if( gate.getName() == source )
                    return gate
            }
            throw new AssertionError("Missing source gate: " + source)
        }

        List<String> feeds(String source) {
            return gate(source).getFeeds().collect { ProcessBacklog it -> it.getName() }.sort()
        }

        /**
         * The times the given source blocked its producer thread
         */
        long blockCount(String source) {
            return gate(source).getBlockCount()
        }

        /**
         * The times the given process paused its operator
         */
        long pauseCount(String name) {
            final backlog = sampler.backlogs.get(name)
            if( backlog == null )
                throw new AssertionError("Missing process backlog: " + name)
            return backlog.getPauseCount()
        }

        WorkflowStats stats() {
            return session.getStatsObserver().getStats()
        }
    }

    /**
     * The last run, also set when it fails
     */
    Run last

    private Session newSession(Map<String,Object> config) {
        final result = new HashMap<String,Object>(config)
        result.put('workDir', workDir.toString())
        return new TestSession(result)
    }

    Run runScript(String text, Map<String,Object> config = [:]) {
        final session = newSession(config)
        session.setBinding(new ScriptBinding())
        session.init(null, null, null, null)
        session.start()

        final loader = ScriptLoaderFactory.create(session)
        // parse(String) is not part of the ScriptLoader interface
        loader.parse(text)
        loader.runScript()
        final result = ScriptHelper.normalizeResult(loader.getResult())

        final sampler = new Sampler(session.getBackpressure())
        // a process registers when it is invoked, and deregisters when it terminates
        for( ProcessBacklog it : session.getBackpressure().processes() )
            sampler.backlogs.put(it.getName(), it)
        last = new Run(result, session, sampler)
        sampler.start()
        try {
            session.fireDataflowNetwork(false)
            session.await()
        }
        finally {
            sampler.done = true
            sampler.join()
            session.destroy()
            // like with GPars, an aborted run does not wait for its operators to stop nor for the
            // process bodies in flight, which may still create task work dirs: wait for both before
            // checking the backlogs and deleting the work dir
            if( session.isAborted() )
                BackpressureInternals.await { session.getBackpressure().processes().isEmpty()
                    && sampler.backlogs.values().every { ProcessBacklog backlog -> backlog.getPreparing() == 0 } }
        }
        assertBalanced(last)
        if( session.getError() != null )
            throw session.getError()
        return last
    }

    /**
     * Check that no gate is left closed and, for a successful run, that every process backlog
     * is back to zero. An aborted run drops the pending tasks and the unconsumed inputs, and
     * relies on releasing every gate instead.
     */
    static void assertBalanced(Run run) {
        final bp = run.backpressure()
        assert !run.sampler.backlogs.isEmpty()
        assert bp.processes().isEmpty()
        for( Backpressure.Gate gate : bp.gates() )
            assert gate.getWaiters() == 0 : "${gate.getName()} waiters"
        if( run.session.getError() != null ) {
            assert bp.isReleased()
            return
        }
        for( ProcessBacklog it : run.sampler.backlogs.values() ) {
            assert it.getPreparing() == 0 : "${it.getName()} preparing"
            assert it.getPending() == 0 : "${it.getName()} pending"
            assert it.ready(Integer.MAX_VALUE) == 0 : "${it.getName()} ready"
        }
        for( Backpressure.Gate gate : bp.gates() )
            assert gate.isOpen(bp.getHighWater()) : "${gate.getName()} open"
    }

    static int intValue(Object value) {
        return ((Number) value).intValue()
    }

    static List<Integer> intValues(Object value) {
        return ((List) value).collect { intValue(it) }
    }

    static String script(String text) {
        return text.replace('$N', String.valueOf(N))
    }

    def 'should bound a slow consumer behind an operator' () {
        when:
        def run = runScript(script('''
            process SLOW {
              input: val x
              output: val x
              exec: sleep 1
            }

            workflow {
              channel.of(1..$N).map { it * 2 } | SLOW | count
            }
            '''))
        then:
        intValue(run.value()) == N
        run.feeds('channel.of') == ['SLOW']
        run.maxBusy('SLOW') <= MAX_BUSY
        // items in transit through the operator are not counted, so the overshoot depends on timing;
        // the block count below is what proves the source was throttled
        run.maxReady('SLOW') < N
        run.pauseCount('SLOW') > 0
        run.blockCount('channel.of') > 0
        // a terminated process is deregistered
        run.backpressure().processes().isEmpty()
    }

    def 'should bound a slow typed process' () {
        when:
        def run = runScript(script('''
            nextflow.enable.types = true

            process SLOW {
              input:
              x: Integer

              output:
              y: Integer

              exec:
              sleep 1
              y = x
            }

            workflow {
              SLOW(channel.of(1..$N)).count()
            }
            '''))
        then:
        intValue(run.value()) == N
        run.feeds('channel.of') == ['SLOW']
        run.maxBusy('SLOW') <= MAX_BUSY
        run.maxReady('SLOW') < N.intdiv(4)
        run.pauseCount('SLOW') > 0
        run.blockCount('channel.of') > 0

        when:
        // the DAG connects the source to the process
        def edges = run.session.getDag().getEdges()
            .findAll { e -> e.getFrom() != null && 'channel.of' == e.getFrom().getLabel() }
        then:
        edges.size() == 1
        edges[0].getTo().getLabel() == 'SLOW'
    }

    def 'should bound a fan-out with a slow branch' () {
        when:
        def run = runScript(script('''
            process FAST {
              input: val x
              output: val x
              exec: true
            }

            process SLOW {
              input: val x
              output: val x
              exec: sleep 1
            }

            workflow {
              ch = channel.of(1..$N)
              FAST(ch)
              SLOW(ch)
              FAST.out.count().combine(SLOW.out.count())
            }
            '''))
        then:
        intValues(run.value()) == [N, N]
        run.feeds('channel.of') == ['FAST', 'SLOW']
        run.maxBusy('SLOW') <= MAX_BUSY
        run.maxReady('SLOW') < N.intdiv(4)
        run.blockCount('channel.of') > 0
    }

    def 'should not deadlock a process reading a channel and a value reduced from it' () {
        when:
        def run = runScript(script('''
            process P {
              input:
                val x
                val all
              output: val x
              exec: true
            }

            workflow {
              ch = channel.of(1..$N)
              P(ch, ch.collect()) | count
            }
            '''))
        then:
        intValue(run.value()) == N
    }

    def 'should not deadlock an operator merging a channel with a value reduced from it' () {
        when:
        def run = runScript(script('''
            process P {
              input: val x
              output: val x
              exec: sleep 1
            }

            workflow {
              a = channel.of(1..$N)
              a.merge(a.collect()) { x, all -> x } | P | count
            }
            '''))
        then:
        // the value input is re-read for every item
        intValue(run.value()) == N
    }

    def 'should not deadlock combine and join' () {
        when:
        def run = runScript(script('''
            process P {
              input: val x
              output: val x
              exec: sleep 1
            }

            workflow {
              left = channel.of(1..$N).map { [it, 'l'] }
              right = channel.of(1..$N).map { [it, 'r'] }
              joined = left.join(right)
              combined = channel.of(1..20).combine(channel.of('a','b'))
              P(joined.mix(combined)) | count
            }
            '''))
        then:
        intValue(run.value()) == N + 40
        run.maxBusy('P') <= MAX_BUSY
    }

    def 'should not deadlock a process with value and multiple queue inputs' () {
        when:
        def run = runScript(script('''
            process Q {
              input: val x
              output: val x
              exec: sleep 1
            }

            process P {
              input:
                val v
                val x
                val y
              output: val x
              exec: true
            }

            workflow {
              src = channel.of(1..$N)
              P(channel.value('v'), channel.of(1..$N), Q(src)) | count
            }
            '''))
        then:
        intValue(run.value()) == N
        run.maxBusy('Q') <= MAX_BUSY
    }

    def 'should not deadlock a process with each inputs' () {
        when:
        def run = runScript('''
            process P {
              input:
                val x
                each y
              output: val x
              exec: sleep 1
            }

            workflow {
              P(channel.of(1..100), [1, 2]) | count
            }
            ''')
        then:
        intValue(run.value()) == 200
        run.maxBusy('P') <= MAX_BUSY
    }

    def 'should not deadlock branch and multiMap' () {
        when:
        def run = runScript(script('''
            process SLOW {
              input: val x
              output: val x
              exec: sleep 1
            }

            process FAST {
              input: val x
              output: val x
              exec: true
            }

            workflow {
              ch = channel.of(1..$N).branch { even: it % 2 == 0; odd: true }
              mm = ch.even.multiMap { a: it; b: it }
              SLOW(mm.a.mix(ch.odd))
              FAST(mm.b)
              SLOW.out.count().combine(FAST.out.count())
            }
            '''))
        then:
        intValues(run.value()) == [N, N.intdiv(2)]
        run.feeds('channel.of') == ['FAST', 'SLOW']
        run.maxBusy('SLOW') <= MAX_BUSY
        run.pauseCount('SLOW') > 0
    }

    def 'should not deadlock topics' () {
        when:
        def run = runScript(script('''
            process P {
              input: val x
              output: val x, topic: 'nums'
              exec: true
            }

            process Q {
              input: val x
              output: val x
              exec: sleep 1
            }

            workflow {
              P(channel.of(1..$N))
              Q(channel.topic('nums')) | count
            }
            '''))
        then:
        intValue(run.value()) == N
        // a topic is not a gated source, Q only throttles itself
        !run.backpressure().gates().collect { Backpressure.Gate it -> it.getName() }.contains('channel.topic')
        run.maxBusy('Q') <= MAX_BUSY
        run.pauseCount('Q') > 0
    }

    /**
     * A recursion feeds one item at a time, so there is nothing to throttle.
     */
    def 'should not deadlock a recursion' () {
        given:
        NextflowMeta.getInstance().getPreview().setRecursion(true)

        when:
        def run = runScript('''
            process foo {
              input: val x
              output: val y
              exec: y = x + 1
            }

            workflow {
              foo.recurse(0).times(50)
              foo.out.last()
            }
            ''')
        then:
        intValue(run.value()) == 50

        cleanup:
        NextflowMeta.getInstance().getPreview().setRecursion(false)
    }

    def 'should bound a slow consumer of fromPath' () {
        given:
        def dir = Files.createDirectories(workDir.resolve('files'))
        for( int i = 0; i < N; i++ )
            Files.createFile(dir.resolve("file_${i}.txt"))

        when:
        def run = runScript(script('''
            process SLOW {
              input: val x
              output: val x
              exec: sleep 1
            }

            workflow {
              channel.fromPath('$DIR/*.txt') | SLOW | count
            }
            ''').replace('$DIR', dir.toString()))
        then:
        intValue(run.value()) == N
        run.feeds('channel.fromPath') == ['SLOW']
        run.maxBusy('SLOW') <= MAX_BUSY
        run.maxReady('SLOW') < N.intdiv(4)
    }

    def 'should bound a slow consumer of fromFilePairs' () {
        given:
        def dir = Files.createDirectories(workDir.resolve('pairs'))
        for( int i = 0; i < 100; i++ ) {
            Files.createFile(dir.resolve("s${i}_1.fq"))
            Files.createFile(dir.resolve("s${i}_2.fq"))
        }

        when:
        def run = runScript('''
            process SLOW {
              input: tuple val(id), val(files)
              output: val id
              exec: sleep 1
            }

            workflow {
              channel.fromFilePairs('$DIR/*_{1,2}.fq') | SLOW | count
            }
            '''.replace('$DIR', dir.toString()))
        then:
        intValue(run.value()) == 100
        run.feeds('channel.fromFilePairs') == ['SLOW']
    }

    def 'should terminate interval with #op' () {
        when:
        def run = runScript('''
            process SLOW {
              input: val x
              output: val x
              exec: sleep 2
            }

            workflow {
              channel.interval('1ms') { i -> i < 2000 ? i : nextflow.Channel.STOP }.$OP | SLOW | count
            }
            '''.replace('$OP', op))
        then:
        intValue(run.value()) == 50
        run.maxBusy('SLOW') <= MAX_BUSY

        where:
        op << ['take(50)', 'until { it >= 50 }']
    }

    def 'should not hang when a process fails while a source is blocked' () {
        when:
        runScript(script('''
            process SLOW {
              errorStrategy 'terminate'

              input: val x
              output: val x
              exec:
                sleep 1
                if( x == 20 ) throw new IllegalStateException('boom')
            }

            workflow {
              channel.of(1..$N) | SLOW | count
            }
            '''))
        then:
        thrown(ProcessFailedException)
        last.blockCount('channel.of') > 0
        last.backpressure().isReleased()
    }

    def 'should bound retried tasks' () {
        when:
        def run = runScript(script('''
            process SLOW {
              errorStrategy 'retry'
              maxRetries 1

              input: val x
              output: val x
              exec:
                sleep 1
                if( task.attempt == 1 && x % 10 == 0 ) throw new IllegalStateException('boom')
            }

            workflow {
              channel.of(1..$N) | SLOW | count
            }
            '''))
        then:
        intValue(run.value()) == N
        run.stats().getRetriesCount() == N.intdiv(10)
        run.maxBusy('SLOW') <= MAX_BUSY
        run.pauseCount('SLOW') > 0
        run.blockCount('channel.of') > 0
    }

    def 'should bound ignored tasks' () {
        when:
        def run = runScript(script('''
            process SLOW {
              errorStrategy 'ignore'

              input: val x
              output: val x
              exec:
                sleep 1
                if( x % 10 == 0 ) throw new IllegalStateException('boom')
            }

            workflow {
              channel.of(1..$N) | SLOW | count
            }
            '''))
        then:
        intValue(run.value()) == N - N.intdiv(10)
        run.stats().getIgnoredCount() == N.intdiv(10)
        run.maxBusy('SLOW') <= MAX_BUSY
        run.pauseCount('SLOW') > 0
        run.blockCount('channel.of') > 0
    }

    def 'should bound tasks skipped by when' () {
        when:
        def run = runScript(script('''
            process SLOW {
              input: val x
              output: val x
              when: x % 2 == 0
              exec: sleep 1
            }

            workflow {
              channel.of(1..$N) | SLOW | count
            }
            '''))
        then:
        intValue(run.value()) == N.intdiv(2)
        run.maxBusy('SLOW') <= MAX_BUSY
        run.pauseCount('SLOW') > 0
    }

    def 'should bound cached tasks' () {
        given:
        def text = script('''
            process SLOW {
              input: val x
              output: val x
              exec: sleep 1
            }

            workflow {
              channel.of(1..$N) | SLOW | count
            }
            ''')
        def config = [resume: UUID.randomUUID().toString()] as Map<String,Object>
        SysEnv.push(NXF_DATAFLOW_HIGH_WATER: String.valueOf(HW), NXF_CACHE_DIR: workDir.resolve('cache').toString())

        when:
        def first = runScript(text, config)
        then:
        intValue(first.value()) == N

        when:
        reset()
        def run = runScript(text, config)
        then:
        intValue(run.value()) == N
        run.stats().getCachedCount() == N
        run.maxBusy('SLOW') <= MAX_BUSY

        cleanup:
        SysEnv.pop()
    }

    def 'should bound tasks skipped by storeDir' () {
        given:
        def store = Files.createDirectories(workDir.resolve('store'))
        def text = '''
            process SLOW {
              storeDir '$STORE'

              input: val x
              output: path "out_${x}.txt"
              exec:
                sleep 1
                // a native task writes its outputs in the store dir, since they are not staged out
                file("$STORE/out_${x}.txt").text = x
            }

            workflow {
              channel.of(1..100) | SLOW | count
            }
            '''.replace('$STORE', store.toString())

        when:
        def first = runScript(text)
        then:
        intValue(first.value()) == 100

        when:
        reset()
        def run = runScript(text)
        then:
        intValue(run.value()) == 100
        run.stats().getProcesses().size() == 1
        run.stats().getProcesses()[0].getStored() == 100
        run.maxBusy('SLOW') <= MAX_BUSY
    }

    def 'should not bound the backlog with backpressure disabled' () {
        given:
        SysEnv.push(NXF_DATAFLOW_HIGH_WATER: '0')

        expect:
        newSession([:]).getBackpressure() == null

        cleanup:
        SysEnv.pop()
    }

    def 'should fail on an invalid high-water mark' () {
        given:
        SysEnv.push(NXF_DATAFLOW_HIGH_WATER: '-1')

        when:
        newSession([:])
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('NXF_DATAFLOW_HIGH_WATER')

        cleanup:
        SysEnv.pop()
    }
}
