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
 */

package nextflow.extension

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.DataflowEventListener
import nextflow.Session
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowHelperTest extends Specification {


    def setupSpec() {
        new Session()
    }

    def 'should subscribe handlers'() {

        when:
        DataflowHelper.checkSubscribeHandlers( [:] )
        then:
        thrown(IllegalArgumentException)

        when:
        DataflowHelper.checkSubscribeHandlers( [ onNext:{}] )
        then:
        true

        when:
        DataflowHelper.checkSubscribeHandlers( [ onNext:{}, xxx:{}] )
        then:
        thrown(IllegalArgumentException)

        when:
        DataflowHelper.checkSubscribeHandlers( [ xxx:{}] )
        then:
        thrown(IllegalArgumentException)
    }

    @Unroll
    def 'should split entry' () {
        when:
        def pair = DataflowHelper.makeKey(pivot, entry)
        then:
        pair.keys == keys
        pair.values == values

        where:
        pivot           | entry                         | keys          | values
        [0]             | ['A','B','C','D','E','F']     | ['A']         | ['B','C','D','E','F']
        [0,1]           | ['A','B','C','D','E','F']     | ['A','B']     | ['C','D','E','F']
        [0,2]           | ['A','B','C','D','E','F']     | ['A','C']     | ['B','D','E','F']
        [0,1,4]         | ['A','B','C','D','E','F']     | ['A','B','E'] | ['C','D','F']
        [0]             | 'A'                           | ['A']         | []
    }

    def 'should validate reduce params' () {
        given:
        def source = new DataflowQueue()
        def target = new DataflowVariable()
        def action = {-> 1}
        def beforeBind = {-> 2}
        def params = new DataflowHelper.ReduceParams()
            .withSource(source)
            .withTarget(target)
            .withSeed('xyz')
            .withAction(action)
            .withBeforeBind(beforeBind)

        expect:
        params.source.is(source)
        params.target.is(target)
        params.seed.is('xyz')
        params.action.is(action)
        params.beforeBind.is(beforeBind)
    }

    def 'should validate operator params' () {
        when:
        def p1 = new DataflowHelper.OpParams().toMap()
        then:
        p1.inputs == List.of()
        p1.outputs == List.of()
        p1.listeners == List.of()

        when:
        def s1 = new DataflowQueue()
        def t1 = new DataflowQueue()
        def l1 = Mock(DataflowEventListener)
        and:
        def p2 = new DataflowHelper.OpParams()
                .withInput(s1)
                .withOutput(t1)
                .withListener(l1)
                .withAccumulator(true)
        then:
        p2.inputs == List.of(s1)
        p2.outputs == List.of(t1)
        p2.listeners == List.of(l1)
        p2.accumulator
        and:
        p2.toMap().inputs == List.of(s1)
        p2.toMap().outputs == List.of(t1)
        p2.toMap().listeners == List.of(l1)

        when:
        def s2 = new DataflowQueue()
        def t2 = new DataflowQueue()
        def l2 = Mock(DataflowEventListener)
        and:
        def p3 = new DataflowHelper.OpParams()
            .withInputs([s1,s2])
            .withOutputs([t1,t2])
            .withListeners([l1,l2])
            .withAccumulator(false)
        then:
        p3.inputs == List.of(s1,s2)
        p3.outputs == List.of(t1,t2)
        p3.listeners == List.of(l1,l2)
        !p3.accumulator
        and:
        p3.toMap().inputs == List.of(s1,s2)
        p3.toMap().outputs == List.of(t1,t2)
        p3.toMap().listeners == List.of(l1,l2)

    }
}
