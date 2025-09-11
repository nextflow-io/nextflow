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


import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class OpClosureTest extends Specification {

    def 'should invoke target closure' () {
        given:
        def code = { a,b -> a+b }
        def context = new ContextSequential()
        def wrapper = new OpClosure(code, context)

        when:
        def result = wrapper.call(1,2)
        then:
        result == 3
        and:
        context.getOperatorRun() != null
    }

    def 'should instrument a closure'() {
        given:
        def code = { int x, int y -> x+y }
        def v1 = 1
        def v2 = 2

        when:
        def c = new OpClosure(code, new ContextSequential())
        def z = c.call([v1, v2] as Object[])
        then:
        z == 3
    }
}
