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

package nextflow.script.params

import static test.TestParser.*

import test.Dsl2Spec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class EnvOutParamTest extends Dsl2Spec {

    def 'should define env outputs' () {
        setup:
        def text = '''
            process hola {
              output:
              env FOO
              env BAR
              
              /echo command/ 
            }
            
            workflow { hola() }
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

    def 'should define env outputs with quotes' () {
        setup:
        def text = '''
            process hola {
              output:
              env 'FOO'
              env 'BAR'
              
              /echo command/ 
            }
            
            workflow { hola() }
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

    def 'should define optional env outputs' () {
        setup:
        def text = '''
            process hola {
              output:
              env FOO optional false
              env BAR optional true

              /echo command/
            }
            
            workflow { hola() }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        EnvOutParam out0 = process.config.getOutputs().get(0)
        EnvOutParam out1 = process.config.getOutputs().get(1)

        then:
        process.config.getOutputs().size() == 2

        out0.getName() == 'FOO'
        out0.getOptional() == false

        out1.getName() == 'BAR'
        out1.getOptional() == true

    }

    def 'should handle invalid env definition' () {
        given:
        def text = '''
            process hola {
              output:
              env { 0 }
              
              /echo command/ 
            }
            
            workflow { hola() }
            '''

        when:
        def binding = [:]
        parseAndReturnProcess(text, binding)

        then:
        def e = thrown(IllegalArgumentException)
        and:
        e.message.startsWith('Unexpected environment output definition')

    }
}
