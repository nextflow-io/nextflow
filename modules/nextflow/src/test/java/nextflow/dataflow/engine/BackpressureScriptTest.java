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
package nextflow.dataflow.engine;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.TimeUnit;

import groovyx.gpars.dataflow.DataflowReadChannel;
import nextflow.Global;
import nextflow.NF;
import nextflow.NextflowMeta;
import nextflow.Session;
import nextflow.SysEnv;
import nextflow.exception.ProcessFailedException;
import nextflow.executor.Executor;
import nextflow.executor.ExecutorFactory;
import nextflow.executor.local.NativeTaskHandler;
import nextflow.processor.TaskHandler;
import nextflow.processor.TaskMonitor;
import nextflow.processor.TaskPollingMonitor;
import nextflow.processor.TaskProcessor;
import nextflow.processor.TaskRun;
import nextflow.script.ScriptBinding;
import nextflow.script.ScriptLoaderFactory;
import nextflow.script.ScriptMeta;
import nextflow.script.ScriptType;
import nextflow.trace.ProgressRecord;
import nextflow.trace.WorkflowStats;
import nextflow.util.Duration;
import org.codehaus.groovy.runtime.InvokerHelper;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import test.ScriptHelper;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

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
class BackpressureScriptTest {

    static final int HW = 8;

    static final int N = 1000;

    /**
     * A task is counted both as preparing and pending between its scheduling and the end of
     * its preparation, so the sampled count may exceed the mark by the tasks being prepared.
     */
    static final int MAX_BUSY = 2 * HW;

    @TempDir
    Path workDir;

    @BeforeEach
    void setup() {
        reset();
        SysEnv.push(Map.of("NXF_DATAFLOW_HIGH_WATER", String.valueOf(HW)));
    }

    private static void reset() {
        TaskProcessor.reset();
        ScriptMeta.reset();
        Global.reset();
        NF.init();
    }

    @AfterEach
    void cleanup() {
        SysEnv.pop();
    }

    static class TestExecutor extends Executor {
        @Override
        protected TaskMonitor createTaskMonitor() {
            return TaskPollingMonitor.create(getSession(), getConfig(), getName(), 2, Duration.of("5ms"));
        }

        @Override
        public TaskHandler createTaskHandler(TaskRun task) {
            // an anonymous subclass, since the constructor is protected
            return new NativeTaskHandler(task, this) {};
        }
    }

    static class TestExecutorFactory extends ExecutorFactory {
        @Override
        protected Class<? extends Executor> getExecutorClass(String executorName) {
            return TestExecutor.class;
        }

        @Override
        protected boolean isTypeSupported(ScriptType type, Object executor) {
            return true;
        }
    }

    static class TestSession extends Session {
        TestSession(Map<String,Object> config) {
            super(config);
        }

        @Override
        public Session start() {
            setExecutorFactory(new TestExecutorFactory());
            return super.start();
        }
    }

    /**
     * Samples the largest backlog of each process during a run, and keeps each backlog
     * since a process is deregistered when it terminates.
     */
    static class Sampler extends Thread {
        final Backpressure bp;
        final Map<String,ProcessBacklog> backlogs = new ConcurrentHashMap<>();
        final Map<String,Integer> busy = new ConcurrentHashMap<>();
        final Map<String,Integer> ready = new ConcurrentHashMap<>();
        volatile boolean done;

        Sampler(Backpressure bp) {
            this.bp = bp;
            setDaemon(true);
        }

        @Override
        public void run() {
            while( !done ) {
                for( var it : bp.processes() ) {
                    backlogs.putIfAbsent(it.getName(), it);
                    busy.merge(it.getName(), it.getPreparing() + it.getPending(), Math::max);
                    ready.merge(it.getName(), it.ready(100_000), Math::max);
                }
                try {
                    Thread.sleep(1);
                }
                catch( InterruptedException e ) {
                    return;
                }
            }
        }
    }

    record Run(Object result, Session session, Sampler sampler) {

        Backpressure backpressure() {
            return session.getBackpressure();
        }

        Object value() throws InterruptedException {
            return ((DataflowReadChannel<?>) result).getVal();
        }

        int maxBusy(String name) {
            return sampler.busy.getOrDefault(name, 0);
        }

        int maxReady(String name) {
            return sampler.ready.getOrDefault(name, 0);
        }

        Backpressure.Gate gate(String source) {
            for( var gate : backpressure().gates() ) {
                if( gate.getName().equals(source) )
                    return gate;
            }
            throw new AssertionError("Missing source gate: " + source);
        }

        List<String> feeds(String source) {
            return gate(source).getFeeds().stream().map(ProcessBacklog::getName).sorted().toList();
        }

        /**
         * The times the given source blocked its producer thread
         */
        long blockCount(String source) {
            return gate(source).getBlockCount();
        }

        /**
         * The times the given process paused its operator
         */
        long pauseCount(String name) {
            final var backlog = sampler.backlogs.get(name);
            if( backlog == null )
                throw new AssertionError("Missing process backlog: " + name);
            return backlog.getPauseCount();
        }

        WorkflowStats stats() {
            return session.getStatsObserver().getStats();
        }
    }

    /**
     * The last run, also set when it fails
     */
    private Run last;

    private Session newSession(Map<String,Object> config) {
        final var result = new HashMap<String,Object>(config);
        result.put("workDir", workDir.toString());
        return new TestSession(result);
    }

    Run run(String text) throws Throwable {
        return run(text, Map.of());
    }

    Run run(String text, Map<String,Object> config) throws Throwable {
        final var session = newSession(config);
        session.setBinding(new ScriptBinding());
        session.init(null, null, null, null);
        session.start();

        final var loader = ScriptLoaderFactory.create(session);
        // parse(String) is not part of the ScriptLoader interface
        InvokerHelper.invokeMethod(loader, "parse", text);
        loader.runScript();
        final var result = ScriptHelper.normalizeResult(loader.getResult());

        final var sampler = new Sampler(session.getBackpressure());
        // a process registers when it is invoked, and deregisters when it terminates
        for( var it : session.getBackpressure().processes() )
            sampler.backlogs.put(it.getName(), it);
        last = new Run(result, session, sampler);
        sampler.start();
        try {
            session.fireDataflowNetwork(false);
            session.await();
        }
        finally {
            sampler.done = true;
            sampler.join();
            session.destroy();
        }
        assertBalanced(last);
        if( session.getError() != null )
            throw session.getError();
        return last;
    }

    /**
     * Check that no gate is left closed and, for a successful run, that every process backlog
     * is back to zero. An aborted run drops the pending tasks and the unconsumed inputs, and
     * relies on releasing every gate instead.
     */
    static void assertBalanced(Run run) {
        final var bp = run.backpressure();
        assertThat(run.sampler().backlogs).isNotEmpty();
        assertThat(bp.processes()).isEmpty();
        for( var gate : bp.gates() )
            assertThat(gate.getWaiters()).as("%s waiters", gate.getName()).isZero();
        if( run.session().getError() != null ) {
            assertThat(bp.isReleased()).isTrue();
            return;
        }
        for( var it : run.sampler().backlogs.values() ) {
            assertThat(it.getPreparing()).as("%s preparing", it.getName()).isZero();
            assertThat(it.getPending()).as("%s pending", it.getName()).isZero();
            assertThat(it.ready(Integer.MAX_VALUE)).as("%s ready", it.getName()).isZero();
        }
        for( var gate : bp.gates() )
            assertThat(gate.isOpen(bp.getHighWater())).as("%s open", gate.getName()).isTrue();
    }

    static int intValue(Object value) {
        return ((Number) value).intValue();
    }

    static List<Integer> intValues(Object value) {
        return ((List<?>) value).stream().map(BackpressureScriptTest::intValue).toList();
    }

    static String script(String text) {
        return text.replace("$N", String.valueOf(N));
    }

    @Test
    void shouldBoundASlowConsumerBehindAnOperator() throws Throwable {
        var run = run(script("""
            process SLOW {
              input: val x
              output: val x
              exec: sleep 1
            }

            workflow {
              channel.of(1..$N).map { it * 2 } | SLOW | count
            }
            """));
        assertThat(intValue(run.value())).isEqualTo(N);
        assertThat(run.feeds("channel.of")).containsExactly("SLOW");
        assertThat(run.maxBusy("SLOW")).isLessThanOrEqualTo(MAX_BUSY);
        // items in transit through the operator are not counted, so the overshoot depends on timing;
        // the block count below is what proves the source was throttled
        assertThat(run.maxReady("SLOW")).isLessThan(N);
        assertThat(run.pauseCount("SLOW")).isPositive();
        assertThat(run.blockCount("channel.of")).isPositive();
        // a terminated process is deregistered
        assertThat(run.backpressure().processes()).isEmpty();
    }

    @Test
    void shouldBoundASlowTypedProcess() throws Throwable {
        var run = run(script("""
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
            """));
        assertThat(intValue(run.value())).isEqualTo(N);
        assertThat(run.feeds("channel.of")).containsExactly("SLOW");
        assertThat(run.maxBusy("SLOW")).isLessThanOrEqualTo(MAX_BUSY);
        assertThat(run.maxReady("SLOW")).isLessThan(N / 4);
        assertThat(run.pauseCount("SLOW")).isPositive();
        assertThat(run.blockCount("channel.of")).isPositive();
        // the DAG connects the source to the process
        var edges = run.session().getDag().getEdges().stream()
            .filter(e -> e.getFrom() != null && "channel.of".equals(e.getFrom().getLabel()))
            .toList();
        assertThat(edges).hasSize(1);
        assertThat(edges.get(0).getTo().getLabel()).isEqualTo("SLOW");
    }

    @Test
    void shouldBoundAFanOutWithASlowBranch() throws Throwable {
        var run = run(script("""
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
            """));
        assertThat(intValues(run.value())).containsExactly(N, N);
        assertThat(run.feeds("channel.of")).containsExactly("FAST", "SLOW");
        assertThat(run.maxBusy("SLOW")).isLessThanOrEqualTo(MAX_BUSY);
        assertThat(run.maxReady("SLOW")).isLessThan(N / 4);
        assertThat(run.blockCount("channel.of")).isPositive();
    }

    @Test
    void shouldNotDeadlockAProcessReadingAChannelAndAValueReducedFromIt() throws Throwable {
        var run = run(script("""
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
            """));
        assertThat(intValue(run.value())).isEqualTo(N);
    }

    @Test
    void shouldNotDeadlockAnOperatorMergingAChannelWithAValueReducedFromIt() throws Throwable {
        var run = run(script("""
            process P {
              input: val x
              output: val x
              exec: sleep 1
            }

            workflow {
              a = channel.of(1..$N)
              a.merge(a.collect()) { x, all -> x } | P | count
            }
            """));
        // the value input is re-read for every item
        assertThat(intValue(run.value())).isEqualTo(N);
    }

    @Test
    void shouldNotDeadlockCombineAndJoin() throws Throwable {
        var run = run(script("""
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
            """));
        assertThat(intValue(run.value())).isEqualTo(N + 40);
        assertThat(run.maxBusy("P")).isLessThanOrEqualTo(MAX_BUSY);
    }

    @Test
    void shouldNotDeadlockAProcessWithValueAndMultipleQueueInputs() throws Throwable {
        var run = run(script("""
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
            """));
        assertThat(intValue(run.value())).isEqualTo(N);
        assertThat(run.maxBusy("Q")).isLessThanOrEqualTo(MAX_BUSY);
    }

    @Test
    void shouldNotDeadlockAProcessWithEachInputs() throws Throwable {
        var run = run("""
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
            """);
        assertThat(intValue(run.value())).isEqualTo(200);
        assertThat(run.maxBusy("P")).isLessThanOrEqualTo(MAX_BUSY);
    }

    @Test
    void shouldNotDeadlockBranchAndMultiMap() throws Throwable {
        var run = run(script("""
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
            """));
        assertThat(intValues(run.value())).containsExactly(N, N / 2);
        assertThat(run.feeds("channel.of")).containsExactly("FAST", "SLOW");
        assertThat(run.maxBusy("SLOW")).isLessThanOrEqualTo(MAX_BUSY);
        assertThat(run.pauseCount("SLOW")).isPositive();
        assertThat(run.blockCount("channel.of")).isPositive();
    }

    @Test
    void shouldNotDeadlockTopics() throws Throwable {
        var run = run(script("""
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
            """));
        assertThat(intValue(run.value())).isEqualTo(N);
        // a topic is not a gated source, Q only throttles itself
        assertThat(run.backpressure().gates()).extracting(Backpressure.Gate::getName).doesNotContain("channel.topic");
        assertThat(run.maxBusy("Q")).isLessThanOrEqualTo(MAX_BUSY);
        assertThat(run.pauseCount("Q")).isPositive();
    }

    /**
     * A recursion feeds one item at a time, so there is nothing to throttle.
     */
    @Test
    void shouldNotDeadlockARecursion() throws Throwable {
        NextflowMeta.getInstance().getPreview().setRecursion(true);
        try {
            var run = run("""
                process foo {
                  input: val x
                  output: val y
                  exec: y = x + 1
                }

                workflow {
                  foo.recurse(0).times(50)
                  foo.out.last()
                }
                """);
            assertThat(intValue(run.value())).isEqualTo(50);
        }
        finally {
            NextflowMeta.getInstance().getPreview().setRecursion(false);
        }
    }

    @Test
    void shouldBoundASlowConsumerOfFromPath() throws Throwable {
        var dir = Files.createDirectories(workDir.resolve("files"));
        for( int i = 0; i < N; i++ )
            Files.createFile(dir.resolve("file_" + i + ".txt"));

        var run = run(script("""
            process SLOW {
              input: val x
              output: val x
              exec: sleep 1
            }

            workflow {
              channel.fromPath('$DIR/*.txt') | SLOW | count
            }
            """).replace("$DIR", dir.toString()));
        assertThat(intValue(run.value())).isEqualTo(N);
        assertThat(run.feeds("channel.fromPath")).containsExactly("SLOW");
        assertThat(run.maxBusy("SLOW")).isLessThanOrEqualTo(MAX_BUSY);
        assertThat(run.maxReady("SLOW")).isLessThan(N / 4);
    }

    @Test
    void shouldBoundASlowConsumerOfFromFilePairs() throws Throwable {
        var dir = Files.createDirectories(workDir.resolve("pairs"));
        for( int i = 0; i < 100; i++ ) {
            Files.createFile(dir.resolve("s" + i + "_1.fq"));
            Files.createFile(dir.resolve("s" + i + "_2.fq"));
        }

        var run = run("""
            process SLOW {
              input: tuple val(id), val(files)
              output: val id
              exec: sleep 1
            }

            workflow {
              channel.fromFilePairs('$DIR/*_{1,2}.fq') | SLOW | count
            }
            """.replace("$DIR", dir.toString()));
        assertThat(intValue(run.value())).isEqualTo(100);
        assertThat(run.feeds("channel.fromFilePairs")).containsExactly("SLOW");
    }

    @ParameterizedTest
    @ValueSource(strings = { "take(50)", "until { it >= 50 }" })
    void shouldTerminateIntervalWith(String op) throws Throwable {
        var run = run("""
            process SLOW {
              input: val x
              output: val x
              exec: sleep 2
            }

            workflow {
              channel.interval('1ms') { i -> i < 2000 ? i : nextflow.Channel.STOP }.$OP | SLOW | count
            }
            """.replace("$OP", op));
        assertThat(intValue(run.value())).isEqualTo(50);
        assertThat(run.maxBusy("SLOW")).isLessThanOrEqualTo(MAX_BUSY);
    }

    @Test
    void shouldNotHangWhenAProcessFailsWhileASourceIsBlocked() {
        assertThatThrownBy(() -> run(script("""
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
            """))).isInstanceOf(ProcessFailedException.class);
        assertThat(last.blockCount("channel.of")).isPositive();
        assertThat(last.backpressure().isReleased()).isTrue();
    }

    @Test
    void shouldBoundRetriedTasks() throws Throwable {
        var run = run(script("""
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
            """));
        assertThat(intValue(run.value())).isEqualTo(N);
        assertThat(run.stats().getRetriesCount()).isEqualTo(N / 10);
        assertThat(run.maxBusy("SLOW")).isLessThanOrEqualTo(MAX_BUSY);
        assertThat(run.pauseCount("SLOW")).isPositive();
        assertThat(run.blockCount("channel.of")).isPositive();
    }

    @Test
    void shouldBoundIgnoredTasks() throws Throwable {
        var run = run(script("""
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
            """));
        assertThat(intValue(run.value())).isEqualTo(N - N / 10);
        assertThat(run.stats().getIgnoredCount()).isEqualTo(N / 10);
        assertThat(run.maxBusy("SLOW")).isLessThanOrEqualTo(MAX_BUSY);
        assertThat(run.pauseCount("SLOW")).isPositive();
        assertThat(run.blockCount("channel.of")).isPositive();
    }

    @Test
    void shouldBoundTasksSkippedByWhen() throws Throwable {
        var run = run(script("""
            process SLOW {
              input: val x
              output: val x
              when: x % 2 == 0
              exec: sleep 1
            }

            workflow {
              channel.of(1..$N) | SLOW | count
            }
            """));
        assertThat(intValue(run.value())).isEqualTo(N / 2);
        assertThat(run.maxBusy("SLOW")).isLessThanOrEqualTo(MAX_BUSY);
        assertThat(run.pauseCount("SLOW")).isPositive();
    }

    @Test
    void shouldBoundCachedTasks() throws Throwable {
        var text = script("""
            process SLOW {
              input: val x
              output: val x
              exec: sleep 1
            }

            workflow {
              channel.of(1..$N) | SLOW | count
            }
            """);
        var config = Map.<String,Object>of("resume", UUID.randomUUID().toString());
        SysEnv.push(Map.of("NXF_DATAFLOW_HIGH_WATER", String.valueOf(HW), "NXF_CACHE_DIR", workDir.resolve("cache").toString()));
        try {
            assertThat(intValue(run(text, config).value())).isEqualTo(N);

            reset();
            var run = run(text, config);
            assertThat(intValue(run.value())).isEqualTo(N);
            assertThat(run.stats().getCachedCount()).isEqualTo(N);
            assertThat(run.maxBusy("SLOW")).isLessThanOrEqualTo(MAX_BUSY);
        }
        finally {
            SysEnv.pop();
        }
    }

    @Test
    void shouldBoundTasksSkippedByStoreDir() throws Throwable {
        var store = Files.createDirectories(workDir.resolve("store"));
        var text = """
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
            """.replace("$STORE", store.toString());
        assertThat(intValue(run(text).value())).isEqualTo(100);

        reset();
        var run = run(text);
        assertThat(intValue(run.value())).isEqualTo(100);
        assertThat(run.stats().getProcesses()).singleElement().extracting(ProgressRecord::getStored).isEqualTo(100);
        assertThat(run.maxBusy("SLOW")).isLessThanOrEqualTo(MAX_BUSY);
    }

    @Test
    void shouldNotBoundTheBacklogWithBackpressureDisabled() {
        SysEnv.push(Map.of("NXF_DATAFLOW_HIGH_WATER", "0"));
        try {
            assertThat(newSession(Map.of()).getBackpressure()).isNull();
        }
        finally {
            SysEnv.pop();
        }
    }

    @Test
    void shouldFailOnAnInvalidHighWaterMark() {
        SysEnv.push(Map.of("NXF_DATAFLOW_HIGH_WATER", "-1"));
        try {
            assertThatThrownBy(() -> newSession(Map.of()))
                .isInstanceOf(IllegalArgumentException.class)
                .hasMessageContaining("NXF_DATAFLOW_HIGH_WATER");
        }
        finally {
            SysEnv.pop();
        }
    }
}
