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

package nextflow.script.params

import static test.TestParser.parseAndReturnProcess

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class EnvOutParamTest extends Specification {

    def 'should define env outputs' () {
        setup:
        def text = '''
            process hola {
              output:
              env FOO into x 
              env BAR into y 
              
              /echo command/ 
            }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        def outs = process.config.getOutputs() as List<EnvOutParam>

        then:
        outs.size() == 2
        and:
        outs[0].name == 'FOO'
        and:
        outs[1].name == 'BAR'

    }

}
