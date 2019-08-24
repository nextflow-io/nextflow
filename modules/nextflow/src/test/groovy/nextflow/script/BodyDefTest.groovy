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

import spock.lang.Timeout

import static test.TestParser.parseAndReturnProcess

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class BodyDefTest extends Specification {

    def 'should set script type properly' () {

        when:
        def body = new BodyDef({->'echo foo'}, 'echo foo', section)
        then:
        body.type == expected
        body.isShell == shell
        where:
        section     | expected              | shell
        'exec'      | ScriptType.GROOVY     | false
        'script'    | ScriptType.SCRIPTLET  | false
        'shell'     | ScriptType.SCRIPTLET  | true
        'workflow'  | ScriptType.GROOVY     | false
    }


    def 'should thrown illegal argument exception for invalid section' () {

        when:
        new BodyDef({}, 'echo foo', 'bar')
        then:
        thrown(IllegalArgumentException)

    }

    def 'should return empty set'() {
        when:
        def body = new BodyDef({->'echo foo'}, 'echo foo')
        then:
        body.getValNames() == [] as Set
    }

    def 'should return variable names referenced in task body'( ) {

        setup:
        def text = '''

        String x = 1

        @Field
        String = 'Ciao'

        z = 'str'

        process hola {

          /
          println $x + $y + $z
          /
        }
        '''
        when:
        def process = parseAndReturnProcess(text)
        then:
        process.taskBody.valRefs == [
                new TokenValRef('x', 13, 20),
                new TokenValRef('y', 13, 25),
                new TokenValRef('z', 13, 30) ] as Set

        process.taskBody.getValNames() == ['x','y','z'] as Set
    }

    def 'should return property names referenced in task body'() {

        when:
        def runner = new TestScriptRunner( process: [executor:'nope'] )
        def script =
                '''
                class Foo { def foo() { return [x:1] };  }

                alpha = 1
                params.beta = 2
                params.zeta = new Foo()
                params.gamma = [:]
                params.hola = 'Ciao'
                delta = new Foo()
                x = 'alpha'

                process simpleTask  {
                    input:
                    val x from 'hola'

                    """
                    echo ${alpha}
                    echo ${params.beta}
                    echo ${params?.gamma?.omega}
                    echo ${params.zeta.foo()}
                    echo ${params.zeta.foo().x}
                    echo ${delta.foo().x}
                    echo ${params."$x"}
                    """
                }

                '''
        runner.setScript(script).execute()
        then:
        runner.getScriptObj().getTaskProcessor().getTaskBody().getValNames() == ['alpha', 'params.beta', 'params.gamma.omega', 'params.zeta', 'delta', 'params', 'x'] as Set

    }

}
