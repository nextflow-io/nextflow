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

package nextflow.provenance

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.processor.TaskId
import nextflow.processor.TaskRun
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProvTrackerTest extends Specification {

    def 'should normalize null values' () {
        given:
        def prov = new ProvTracker()
        and:
        def t1 = new TaskRun(id: new TaskId(1))
        def t2 = new TaskRun(id: new TaskId(2))
        and:
        def msg1 = ['foo']
        def msg2 = ['foo', new ProvTracker.NullMessage()]

        when:
        def result1 = prov.beforeRun(t1, msg1)
        then:
        result1 == msg1
        and:
        t1.upstreamTasks == [] as Set

        when:
        def result2 = prov.beforeRun(t2, msg2)
        then:
        result2 == ['foo', null]
        and:
        t2.upstreamTasks == [] as Set
    }

    def 'should bind value to task run' () {
        given:
        def prov = new ProvTracker()
        and:
        def t1 = new TaskRun(id: new TaskId(1))
        def c1 = new DataflowQueue()
        def v1 = 'foo'
        
        when:
        prov.bindOutput(t1, c1, v1)
        then:
        c1.val == 'foo'
        and:
        prov.@messages.get(System.identityHashCode(v1)) == t1
    }

    def 'should determine upstream tasks' () {
        given:
        def prov = new ProvTracker()
        and:
        def t1 = new TaskRun(id: new TaskId(1))
        def t2 = new TaskRun(id: new TaskId(2))
        def t3 = new TaskRun(id: new TaskId(3))
        and:
        def v1 = new Object()
        def v2 = new Object()
        def v3 = new Object()

        when:
        prov.bindOutput(t1, Mock(DataflowWriteChannel), v1)
        prov.bindOutput(t2, Mock(DataflowWriteChannel), v2)
        and:
        prov.beforeRun(t3, [v1, v2])
        then:
        t3.upstreamTasks == [t1.id, t2.id] as Set
    }

}
