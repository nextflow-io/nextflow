/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskClosureTest extends Specification {

    def 'should invoke inner closure' () {

        given:
        def code = { String it -> it.reverse() }

        when:
        def wrap = new TaskClosure(code, 'reverse')
        then:
        wrap.getMaximumNumberOfParameters() == 1
        wrap.getParameterTypes() as List == [String]
        wrap.call('Hello') == 'olleH'

    }

    def 'should clone and invoke the closure' ()  {

        given:
        def code = { String a, Integer b -> a.toInteger() + b }
        def wrap = new TaskClosure(code, 'closure source code')

        when:
        def copy = (TaskClosure)wrap.clone()
        then:
        copy.getSource() == 'closure source code'
        !copy.getInnerClosure().is( code )
        copy.getMaximumNumberOfParameters() == 2
        copy.getParameterTypes() as List == [String, Integer]
        copy.call('3', 2) == 5

    }
}
