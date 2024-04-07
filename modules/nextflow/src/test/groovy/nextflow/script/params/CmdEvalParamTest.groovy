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
class CmdEvalParamTest extends Dsl2Spec {

    def 'should define eval outputs' () {
        setup:
        def text = '''
            process hola {
              output:
              eval 'foo --version' 
              eval "$params.cmd --help"
              eval "$tool --test"  
              
              /echo command/ 
            }
            
            workflow { hola() }
            '''

        def binding = [params:[cmd:'bar'], tool: 'other']
        def process = parseAndReturnProcess(text, binding)

        when:
        def outs = process.config.getOutputs() as List<CmdEvalParam>

        then:
        outs.size() == 3
        and:
        outs[0].getName() =~ /nxf_out_eval_\d+/
        outs[0].getTarget(binding) == 'foo --version'
        and:
        outs[1].getName() =~ /nxf_out_eval_\d+/
        outs[1].getTarget(binding) == 'bar --help'
        and:
        outs[2].getName() =~ /nxf_out_eval_\d+/
        outs[2].getTarget(binding) == 'other --test'
    }

}
