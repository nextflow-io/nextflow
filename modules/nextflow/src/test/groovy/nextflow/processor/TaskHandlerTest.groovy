/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
import nextflow.executor.Executor
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
                memory: '4 GB'
        ]
        def task = new TaskRun(id: new TaskId(100), workDir: folder, name:'task1', exitStatus: 127, config: config  )
        task.metaClass.getHashLog = { "5d5d7ds" }
        task.processor = Mock(TaskProcessor)
        task.processor.getSession() >> new Session()
        task.processor.getName() >> 'TheProcessName'
        task.processor.getExecutor() >> Mock(Executor)
        task.processor.getProcessEnvironment() >> [FOO:'hola', BAR: 'mundo', AWS_SECRET: '12345']
        task.context = new TaskContext(Mock(Script), [:], 'none')

        def handler = [:] as TaskHandler
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
        handler = [:] as TaskHandler
        handler.status = TaskStatus.COMPLETED
        handler.submitTimeMillis = 1000
        handler.startTimeMillis = 1500
        handler.task = task.clone()
        handler.task.failed = true
        then:
        handler.getTraceRecord().status == 'FAILED'


        when:
        handler = [:] as TaskHandler
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

}
