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

package nextflow.script

import spock.lang.Specification

import nextflow.extension.OpCall

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WorkflowBindingTest extends Specification {

    def setupSpec() {
        WorkflowBinding.init()
    }

    def 'should lookup variables' () {
        given:
        def obj1 = new Object()
        def obj2 = new Object()
        def binding = new WorkflowBinding()

        when:
        binding.foo = obj1
        binding.bar = obj2
        then:
        WorkflowBinding.lookup(obj1) == 'foo'
        WorkflowBinding.lookup(obj2) == 'bar'
    }

    def 'should invoke an extension method' () {

        given:
        def FOO = Mock(ComponentDef)
        def ARGS = ['alpha','beta'] as Object[]
        def binding = Spy(WorkflowBinding)
        binding.@meta = Mock(ScriptMeta)
        
        // should invoke an extension component
        when:
        def result = binding.invokeMethod('foo', ARGS)
        then:
        1 * binding.getComponent0('foo') >> FOO
        1 * FOO.invoke_o(ARGS) >> 'Hello'
        result == 'Hello'

        // should invoke an extension operator
        when:
        result = binding.invokeMethod('map', ARGS)
        then:
        1 * binding.getComponent0('map') >> null
        result instanceof OpCall
        (result as OpCall).methodName == 'map'
        (result as OpCall).args == ARGS

        // should throw missing method exception 
        when:
        binding.invokeMethod('foo', ARGS)
        then:
        1 * binding.getComponent0('foo') >> null
        thrown(MissingMethodException)

    }

    def 'should get an extension method as variable' () {

        given:
        def FOO = Mock(ComponentDef)
        def binding = Spy(WorkflowBinding)
        binding.setVariable('alpha', 'Hello')
        binding.@meta = Mock(ScriptMeta)

        // should return an existing variable
        when:
        def result = binding.getVariable('alpha')
        then:
        0 * binding.getComponent0('alpha') >> null
        result == 'Hello'

        when:
        result = binding.getVariable('foo')
        then:
        1 * binding.getComponent0('foo') >> FOO
        result == FOO

        when:
        result = binding.getVariable('map')
        then:
        1 * binding.getComponent0('map') >> null
        result instanceof OpCall
        (result as OpCall).methodName == 'map'
        (result as OpCall).args == [] as Object[]

    }

}
