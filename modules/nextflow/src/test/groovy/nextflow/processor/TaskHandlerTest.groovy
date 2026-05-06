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

package nextflow.processor

import java.util.concurrent.atomic.LongAdder

import nextflow.Session
import nextflow.SysEnv
import nextflow.executor.Executor
import nextflow.trace.TraceRecord
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.lang.Unroll
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskHandlerTest extends Specification {

    def static final long KB = 1024

    def 'test get trace record'() {

        given:
        def traceText =  '''
        pid state %cpu %mem vmem rss peak_vmem peak_rss rchar wchar syscr syscw read_bytes write_bytes
        1 0 10 20 11084 1220 21084 2220 4790 12 11 1 20 30
        '''
                .leftTrim()

        def folder = TestHelper.createInMemTempDir()
        folder.resolve( TaskRun.CMD_TRACE ).text = traceText

        def config = [
                tag: 'seq_x',
                container: 'ubuntu',
                queue: 'longjobs',
                cpus: 2,
                time: '1 hour',
                disk: '100 GB',
                memory: '4 GB',
                accelerator: [request: 3, type: 'v100']
        ]
        def task = new TaskRun(id: new TaskId(100), workDir: folder, name:'task1', exitStatus: 127, config: config  )
        task.metaClass.getHashLog = { "5d5d7ds" }
        task.processor = Mock(TaskProcessor)
        task.processor.getSession() >> new Session()
        task.processor.getName() >> 'TheProcessName'
        task.processor.getExecutor() >> Mock(Executor)
        task.processor.getProcessEnvironment() >> [FOO:'hola', BAR: 'mundo', AWS_SECRET: '12345']
        task.context = new TaskContext(Mock(Script), [:], 'none')

        def handler = Spy(TaskHandler)
        handler.task = task
        handler.status = TaskStatus.COMPLETED
        handler.submitTimeMillis = 1000
        handler.startTimeMillis = 1500

        when:
        def trace = handler.getTraceRecord()

        then:
        trace.task_id == 100
        trace.status == 'COMPLETED'
        trace.hash == '5d5d7ds'
        trace.name == 'task1'
        trace.exit == 127
        trace.submit == 1000
        trace.start == 1500
        trace.process == 'TheProcessName'
        trace.tag == 'seq_x'
        trace.'%cpu' == 1.0f
        trace.'%mem' == 2.0f
        trace.rss == 1220 * KB
        trace.vmem == 11084 * KB
        trace.peak_rss == 2220 * KB
        trace.peak_vmem == 21084 * KB
        trace.rchar == 4790
        trace.wchar == 12
        trace.syscr == 11
        trace.syscw ==  1
        trace.read_bytes == 20
        trace.write_bytes == 30
        trace.queue == 'longjobs'
        trace.cpus == 2
        trace.time == Duration.of('1 hour').toMillis()
        trace.memory == MemoryUnit.of('4 GB').toBytes()
        trace.disk == MemoryUnit.of('100 GB').toBytes()
        trace.env == 'FOO=hola\nBAR=mundo\nAWS_SECRET=[secure]\n'
        trace.accelerator == 3
        trace.accelerator_type == 'v100'

        // check get method
        trace.getFmtStr('%cpu') == '1.0%'
        trace.getFmtStr('%mem') == '2.0%'
        trace.getFmtStr('rss') == '1.2 MB'
        trace.getFmtStr('vmem') == '10.8 MB'
        trace.getFmtStr('peak_rss') == '2.2 MB'
        trace.getFmtStr('peak_vmem') == '20.6 MB'
        trace.getFmtStr('time') == '1h'
        trace.getFmtStr('memory') == '4 GB'
        trace.getFmtStr('disk') == '100 GB'

        when:
        handler = Spy(TaskHandler)
        handler.status = TaskStatus.COMPLETED
        handler.submitTimeMillis = 1000
        handler.startTimeMillis = 1500
        handler.task = task.clone()
        handler.task.failed = true
        then:
        handler.getTraceRecord().status == 'FAILED'


        when:
        handler = Spy(TaskHandler)
        handler.status = TaskStatus.COMPLETED
        handler.submitTimeMillis = 1000
        handler.startTimeMillis = 1500
        handler.task = task.clone()
        handler.task.aborted = true
        then:
        handler.getTraceRecord().status == 'ABORTED'

    }

    LongAdder _adder(Integer x) {
        if( x != null ) {
            def adder = new LongAdder()
            adder.add(x)
            return adder
        }
        else
            return null
    }

    @Unroll
    def 'should validate has forks' () {
        given:
        def handler = Spy(TaskHandler)
        handler.task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) {
                getMaxForks() >> MAX
                getForksCount() >> { COUNT ? _adder(COUNT) : null }
            }
        }

        expect:
        handler.canForkProcess() == EXPECT

        where:
        COUNT   | MAX   | EXPECT
        null    | null  | true      // when max is null or 0 result is `true by default
        0       | 0     | true
        1       | 0     | true
        and:
        1       | 10    | true
        9       | 10    | true
        10      | 10    | false
        90      | 10    | false
        1       | 1     | false

    }

    def 'should increment forked count' () {
        given:
        def COUNTER = 10
        def adder = new LongAdder(); adder.add(COUNTER)
        def handler = Spy(TaskHandler)
        handler.task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) {
                getForksCount() >> { adder }
            }
        }

        when:
        handler.incProcessForks()
        then:
        handler.task.processor.getForksCount().intValue() == COUNTER +1
    }

    def 'should decrement forked count' () {
        given:
        def COUNTER = 10
        def adder = new LongAdder(); adder.add(COUNTER)
        def handler = Spy(TaskHandler)
        handler.task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) {
                getForksCount() >> { adder }
            }
        }

        when:
        handler.decProcessForks()
        then:
        handler.task.processor.getForksCount().intValue() == COUNTER -1
    }

    def 'should validate is submit timeout' () {
        given:
        def handler = Spy(TaskHandler)
        handler.status = TaskStatus.SUBMITTED
        handler.task = Mock(TaskRun) {
            getConfig() >> Mock(TaskConfig) { getMaxSubmitAwait() >> Duration.of('500ms') }
        }

        when:
        def timeout = handler.isSubmitTimeout()
        then:
        !timeout

        when:
        sleep 1_000
        and:
        timeout = handler.isSubmitTimeout()
        then:
        timeout

    }

    @Unroll
    def 'should validate status' () {
        given:
        def handler = Spy(TaskHandler)
        handler.status = STATUS

        expect:
        handler.isNew() == EXPECT_NEW
        handler.isRunning() == EXPECT_RUNNING
        handler.isSubmitted() == EXPECT_SUBMITTED
        handler.isActive() == EXPECTED_ACTIVE
        handler.isCompleted() == EXPECT_COMPLETE

        where:
        STATUS              | EXPECT_NEW  | EXPECT_SUBMITTED | EXPECT_RUNNING | EXPECTED_ACTIVE | EXPECT_COMPLETE
        TaskStatus.NEW      | true        | false            | false          | false           | false
        TaskStatus.SUBMITTED| false       | true             | false          | true            | false
        TaskStatus.RUNNING  | false       | false            | true           | true            | false
        TaskStatus.COMPLETED| false       | false            | false          | false           | true
    }

    @Unroll
    def 'should include the tower prefix'() {
        given:
        def name = 'job_1'

        expect:
        TaskHandler.prependWorkflowPrefix(name, ENV) == EXPECTED

        where:
        ENV                         | EXPECTED
        [:]                         | "job_1"
        [TOWER_WORKFLOW_ID: '1234'] | "tw-1234-job_1"
    }

    def 'should not kill task twice'() {
        given:
        def handler = Spy(TaskHandler)
        when:
        handler.kill()
        then:
        1 * handler.killTask() >> {}

        when:
        handler.kill()
        then:
        0 * handler.killTask()
    }

    def 'should parse Fusion trace and set gpuMetrics'() {
        given:
        def folder = TestHelper.createInMemTempDir()
        folder.resolve('.fusion').mkdir()
        folder.resolve(TaskRun.FUSION_TRACE).text = '{"proc":{"realtime":100},"gpu":{"name":"Tesla T4","pct":75,"peak":100}}'

        def task = new TaskRun(workDir: folder)
        task.processor = Mock(TaskProcessor)
        task.processor.getExecutor() >> Mock(Executor) { isFusionEnabled() >> true }

        def handler = Spy(TaskHandler)
        handler.task = task

        def record = new TraceRecord()

        when:
        handler.parseFusionTrace(record)

        then:
        record.gpuMetrics.name == 'Tesla T4'
        record.gpuMetrics.pct == 75
        record.gpuMetrics.peak == 100
    }

    def 'should not read Fusion trace when Fusion is disabled'() {
        given:
        def folder = TestHelper.createInMemTempDir()
        def task = new TaskRun(workDir: folder)
        task.processor = Mock(TaskProcessor)
        task.processor.getExecutor() >> Mock(Executor) { isFusionEnabled() >> false }

        def handler = Spy(TaskHandler)
        handler.task = task

        def record = new TraceRecord()

        when:
        handler.parseFusionTrace(record)

        then:
        record.gpuMetrics == null
    }

    def 'should handle missing Fusion trace file gracefully'() {
        given:
        def folder = TestHelper.createInMemTempDir()
        def task = new TaskRun(workDir: folder)
        task.processor = Mock(TaskProcessor)
        task.processor.getExecutor() >> Mock(Executor) { isFusionEnabled() >> true }

        def handler = Spy(TaskHandler)
        handler.task = task

        def record = new TraceRecord()

        when:
        handler.parseFusionTrace(record)

        then:
        noExceptionThrown()
        record.gpuMetrics == null
    }

    def 'should parse full Fusion trace metrics when NXF_FUSION_TRACE is enabled'() {
        given:
        SysEnv.push([NXF_FUSION_TRACE: 'true'])
        def folder = TestHelper.createInMemTempDir()
        folder.resolve('.fusion').mkdir()
        folder.resolve(TaskRun.FUSION_TRACE).text = '{"proc":{"realtime":660541,"pct_cpu":1045,"cpu_name":"Intel CPU","rchar":100,"wchar":200,"syscr":10,"syscw":20,"read_bytes":300,"write_bytes":400,"pct_mem":56,"vmem":1000,"rss":500,"peak_vmem":1100,"peak_rss":600,"vol_ctxt":100,"inv_ctxt":50},"gpu":{"name":"Tesla T4","pct":75,"peak":100},"cgroup":{"version":"v2","memory_current":25469927424,"memory_peak":41178980352,"memory_rss":67919872,"memory_peak_rss":14783070208,"memory_limit":77309411328}}'

        def task = new TaskRun(workDir: folder)
        task.processor = Mock(TaskProcessor)
        task.processor.getExecutor() >> Mock(Executor) { isFusionEnabled() >> true }

        def handler = Spy(TaskHandler)
        handler.task = task

        def record = new TraceRecord()

        when:
        handler.parseFusionTrace(record)

        then:
        // cgroup memory metrics
        record.store.get('vmem') == 25469927424L
        record.store.get('rss') == 67919872L
        record.store.get('peak_vmem') == 41178980352L
        record.store.get('peak_rss') == 14783070208L
        // proc metrics
        record.store.get('realtime') == 660541L
        record.store.get('%cpu') == 104.5f
        record.store.get('cpu_model') == 'Intel CPU'
        // gpu metrics
        record.gpuMetrics.name == 'Tesla T4'
        record.gpuMetrics.pct == 75

        cleanup:
        SysEnv.pop()
    }

    def 'should only parse GPU metrics when NXF_FUSION_TRACE is disabled'() {
        given:
        SysEnv.push([NXF_FUSION_TRACE: 'false'])
        def folder = TestHelper.createInMemTempDir()
        folder.resolve('.fusion').mkdir()
        folder.resolve(TaskRun.FUSION_TRACE).text = '{"proc":{"realtime":660541,"pct_cpu":1045},"gpu":{"name":"Tesla T4","pct":75,"peak":100},"cgroup":{"version":"v2","memory_current":25469927424}}'

        def task = new TaskRun(workDir: folder)
        task.processor = Mock(TaskProcessor)
        task.processor.getExecutor() >> Mock(Executor) { isFusionEnabled() >> true }

        def handler = Spy(TaskHandler)
        handler.task = task

        def record = new TraceRecord()

        when:
        handler.parseFusionTrace(record)

        then:
        // only GPU metrics populated (existing behavior)
        record.gpuMetrics.name == 'Tesla T4'
        // no proc/cgroup metrics populated
        record.store.get('realtime') == null
        record.store.get('vmem') == null

        cleanup:
        SysEnv.pop()
    }

    def 'should skip command trace file when NXF_FUSION_TRACE is enabled'() {
        given:
        SysEnv.push([NXF_FUSION_TRACE: 'true'])
        def folder = TestHelper.createInMemTempDir()
        // Create a .command.trace with different values to prove it's NOT read
        folder.resolve(TaskRun.CMD_TRACE).text = '''\
            nextflow.trace/v2
            realtime=99999
            %cpu=100
            '''.stripIndent().leftTrim()
        // Create fusion trace
        folder.resolve('.fusion').mkdir()
        folder.resolve(TaskRun.FUSION_TRACE).text = '{"proc":{"realtime":1000,"pct_cpu":500,"cpu_name":"CPU","rchar":0,"wchar":0,"syscr":0,"syscw":0,"read_bytes":0,"write_bytes":0,"pct_mem":0,"vmem":0,"rss":0,"peak_vmem":0,"peak_rss":0,"vol_ctxt":0,"inv_ctxt":0}}'

        def task = new TaskRun(workDir: folder)
        task.processor = Mock(TaskProcessor)
        task.processor.getSession() >> new Session()
        task.processor.getExecutor() >> Mock(Executor) { isFusionEnabled() >> true }
        task.processor.getName() >> 'test'
        task.processor.getProcessEnvironment() >> [:]
        task.config = new TaskConfig()
        task.context = new TaskContext(Mock(Script), [:], 'none')

        def handler = Spy(TaskHandler)
        handler.task = task
        handler.status = TaskStatus.COMPLETED
        handler.submitTimeMillis = 500L
        handler.startTimeMillis = 1000L
        handler.completeTimeMillis = 2000L

        when:
        def record = handler.getTraceRecord()

        then:
        // Should use Fusion realtime (1000), NOT .command.trace realtime (99999)
        record.store.get('realtime') == 1000L

        cleanup:
        SysEnv.pop()
    }

    @Unroll
    def 'should set isChildArray flag'() {
        given:
        def handler = Spy(TaskHandler)

        expect:
        !handler.isArrayChild
        and:
        handler.withArrayChild(VALUE).isArrayChild == VALUE

        where:
        VALUE   | _
        false   | _
        true    | _
    }
}
