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

package nextflow.extension.op

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.operator.DataflowEventListener
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class OpTest extends Specification {

    def 'should validate operator params' () {
        when:
        def p1 = new Op().toMap()
        then:
        p1.inputs == List.of()
        p1.outputs == List.of()
        p1.listeners == List.of()

        when:
        def s1 = new DataflowQueue()
        def t1 = new DataflowQueue()
        def l1 = Mock(DataflowEventListener)
        def c1 = { return 'foo' }
        and:
        def p2 = new Op()
            .withInput(s1)
            .withOutput(t1)
            .withListener(l1)
            .withCode(c1)
        then:
        p2.inputs == List.of(s1)
        p2.outputs == List.of(t1)
        p2.listeners == List.of(l1)
        p2.context instanceof ContextSequential
        p2.code == c1
        and:
        p2.toMap().inputs == List.of(s1)
        p2.toMap().outputs == List.of(t1)
        p2.toMap().listeners == List.of(l1)

        when:
        def s2 = new DataflowQueue()
        def t2 = new DataflowQueue()
        def l2 = Mock(DataflowEventListener)
        and:
        def p3 = new Op()
            .withInputs([s1,s2])
            .withOutputs([t1,t2])
            .withListeners([l1,l2])
            .withContext(new ContextGrouping())
        then:
        p3.inputs == List.of(s1,s2)
        p3.outputs == List.of(t1,t2)
        p3.listeners == List.of(l1,l2)
        p3.context instanceof ContextGrouping
        and:
        p3.toMap().inputs == List.of(s1,s2)
        p3.toMap().outputs == List.of(t1,t2)
        p3.toMap().listeners == List.of(l1,l2)
    }

}
