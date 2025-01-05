/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.prov

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.processor.TaskConfig
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TrackerTest extends Specification {

    def 'should normalize null values' () {
        given:
        def prov = new Tracker()
        and:
        def t1 = new TaskRun(id: new TaskId(1), processor: Mock(TaskProcessor), config: Mock(TaskConfig))
        def t2 = new TaskRun(id: new TaskId(2), processor: Mock(TaskProcessor), config: Mock(TaskConfig))
        and:
        def msg1 = [Tracker.Msg.of('foo')]
        def msg2 = [Tracker.Msg.of('foo'), Tracker.Msg.of(null)]

        when:
        def result1 = prov.receiveInputs(t1, msg1)
        then:
        result1 == msg1.value
        and:
        t1.upstreamTasks == [] as Set

        when:
        def result2 = prov.receiveInputs(t2, msg2)
        then:
        result2 == ['foo', null]
        and:
        t2.upstreamTasks == [] as Set
    }

    def 'should bind value to task run' () {
        given:
        def prov = new Tracker()
        and:
        def t1 = new TaskRun(id: new TaskId(1), processor: Mock(TaskProcessor), config: Mock(TaskConfig))
        def c1 = new DataflowQueue()
        def v1 = 'foo'
        
        when:
        def m1 = prov.bindOutput(t1, c1, v1)
        then:
        c1.val.is(m1)
        and:
        prov.@messages.get(System.identityHashCode(m1)) == t1
    }

    def 'should determine upstream tasks' () {
        given:
        def prov = new Tracker()
        and:
        def t1 = new TaskRun(id: new TaskId(1), processor: Mock(TaskProcessor), config: Mock(TaskConfig))
        def t2 = new TaskRun(id: new TaskId(2), processor: Mock(TaskProcessor), config: Mock(TaskConfig))
        def t3 = new TaskRun(id: new TaskId(3), processor: Mock(TaskProcessor), config: Mock(TaskConfig))
        and:
        def v1 = new Object()
        def v2 = new Object()

        when:
        def m1 = prov.bindOutput(t1, Mock(DataflowWriteChannel), v1)
        def m2 = prov.bindOutput(t2, Mock(DataflowWriteChannel), v2)
        and:
        prov.receiveInputs(t3, [m1, m2])
        then:
        t3.upstreamTasks == [t1.id, t2.id] as Set
    }

    def 'should determine upstream task with operator' () {
        given:
        def prov = new Tracker()
        and:
        def v1 = Integer.valueOf(1)
        def v2 = Integer.valueOf(2)
        def v3 = Integer.valueOf(3)
        and:
        def t1 = new TaskRun(id: new TaskId(1), processor: Mock(TaskProcessor), config: Mock(TaskConfig))
        def p2 = new OperatorRun()
        def t3 = new TaskRun(id: new TaskId(3), processor: Mock(TaskProcessor), config: Mock(TaskConfig))

        when:
        prov.receiveInputs(t1, [])
        def m1 = prov.bindOutput(t1, Mock(DataflowWriteChannel), v1)
        and:
        prov.receiveInputs(p2, [m1])
        and:
        def m2 = prov.bindOutput(p2, Mock(DataflowWriteChannel), v2)
        and:
        prov.receiveInputs(t3, [m2])
        and:
        def m3 = prov.bindOutput(t3, Mock(DataflowWriteChannel), v3)

        then:
        t3.upstreamTasks == [t1.id] as Set
    }

}
