/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.script
import static test.TestParser.parseAndReturnProcess

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskBodyTest extends Specification {

    def 'should set script type properly' () {

        when:
        def body = new TaskBody({->'echo foo'}, 'echo foo', section)
        then:
        body.type == expected
        body.isShell == shell
        where:
        section     | expected              | shell
        'exec'      | ScriptType.GROOVY     | false
        'script'    | ScriptType.SCRIPTLET  | false
        'shell'     | ScriptType.SCRIPTLET  | true

    }


    def 'should thrown illegal argument exception for invalid section' () {

        when:
        new TaskBody({}, 'echo foo', 'bar')
        then:
        thrown(IllegalArgumentException)

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
        def runner = new ScriptRunner( process: [executor:'nope'] )
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
