/*
 * Copyright 2013-2023, Seqera Labs
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

import static test.TestParser.*

import test.Dsl2Spec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdOutParamTest extends Dsl2Spec {

    def 'should define env outputs' () {
        setup:
        def text = '''
            process hola {
              output:
              cmd 'foo --version' 
              cmd 'bar --help'
              
              /echo command/ 
            }
            
            workflow { hola() }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        def outs = process.config.getOutputs() as List<CmdOutParam>

        then:
        outs.size() == 2
        and:
        outs[0].name == 'nxf_out_cmd_1'
        outs[0].target == 'foo --version'
        and:
        outs[1].name == 'nxf_out_cmd_2'
        outs[1].target == 'bar --help'

    }

}
