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

import groovyx.gpars.dataflow.DataflowVariable
import test.Dsl2Spec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TupleOutParamTest extends Dsl2Spec {

    def 'should define output tuples'() {

        setup:
        def text = '''
            process hola {
              output:
                tuple val(x)
                tuple val(y), stdout, file('*.fa') 
                tuple stdout, val(z)

              return ''
            }
            
            workflow {
              hola()
            }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        TupleOutParam out0 = process.config.getOutputs().get(0)
        TupleOutParam out1 = process.config.getOutputs().get(1)
        TupleOutParam out2 = process.config.getOutputs().get(2)

        then:
        process.config.getOutputs().size() == 3

        out0.outChannel instanceof DataflowVariable
        out0.inner.size() == 1
        out0.inner[0] instanceof ValueOutParam
        out0.inner[0].name == 'x'
        out0.inner[0].index == 0

        out1.outChannel instanceof DataflowVariable
        out1.inner[0] instanceof ValueOutParam
        out1.inner[0].name == 'y'
        out1.inner[0].index == 1
        out1.inner[1] instanceof StdOutParam
        out1.inner[1].name == '-'
        out1.inner[1].index == 1
        out1.inner[2] instanceof FileOutParam
        out1.inner[2].name == null
        out1.inner[2].filePattern == '*.fa'
        out1.inner[2].index == 1
        out1.inner.size() ==3

        out2.outChannel instanceof DataflowVariable
        out2.inner.size() == 2
        out2.inner[0] instanceof StdOutParam
        out2.inner[0].name == '-'
        out2.inner[0].index == 2
        out2.inner[1] instanceof ValueOutParam
        out2.inner[1].name == 'z'
        out2.inner[1].index == 2

    }

    def 'should define output tuples 2'() {

        setup:
        def text = '''
            process hola {
              output:
                tuple val(x)
                tuple val(y), stdout, file('*.fa')
                tuple stdout, val(z)

              return ''
            }
            
            workflow {
              hola()
            }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        TupleOutParam out0 = process.config.getOutputs().get(0)
        TupleOutParam out1 = process.config.getOutputs().get(1)
        TupleOutParam out2 = process.config.getOutputs().get(2)

        then:
        process.config.getOutputs().size() == 3

        out0.outChannel instanceof DataflowVariable
        out0.inner.size() == 1
        out0.inner[0] instanceof ValueOutParam
        out0.inner[0].name == 'x'
        out0.inner[0].index == 0

        out1.outChannel instanceof DataflowVariable
        out1.inner[0] instanceof ValueOutParam
        out1.inner[0].name == 'y'
        out1.inner[0].index == 1
        out1.inner[1] instanceof StdOutParam
        out1.inner[1].name == '-'
        out1.inner[1].index == 1
        out1.inner[2] instanceof FileOutParam
        out1.inner[2].name == null
        out1.inner[2].filePattern == '*.fa'
        out1.inner[2].index == 1
        out1.inner.size() ==3

        out2.outChannel instanceof DataflowVariable
        out2.inner.size() == 2
        out2.inner[0] instanceof StdOutParam
        out2.inner[0].name == '-'
        out2.inner[0].index == 2
        out2.inner[1] instanceof ValueOutParam
        out2.inner[1].name == 'z'
        out2.inner[1].index == 2

    }

    def 'should create tuple of env' () {
        setup:
        def text = '''
            process hola {
              output:
                tuple env(FOO), env(BAR)
              
              /echo command/ 
            }
            
            workflow {
              hola()
            }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        def outs = process.config.getOutputs() as List<TupleOutParam>
        then:
        println outs.outChannel
        outs.size() == 1
        and:
        outs[0].outChannel instanceof DataflowVariable
        and:
        outs[0].inner.size() ==2
        and:
        outs[0].inner[0] instanceof EnvOutParam
        outs[0].inner[0].getName() == 'FOO'
        and:
        outs[0].inner[1] instanceof EnvOutParam
        outs[0].inner[1].getName() == 'BAR'
    }

    def 'should create tuple of eval' () {
        setup:
        def text = '''
            process hola {
              output:
                tuple eval('this --one'), eval("$other --two")
              
              /echo command/ 
            }
            
            workflow {
              hola()
            }
            '''

        def binding = [other:'tool']
        def process = parseAndReturnProcess(text, binding)

        when:
        def outs = process.config.getOutputs() as List<TupleOutParam>
        then:
        println outs.outChannel
        outs.size() == 1
        and:
        outs[0].outChannel instanceof DataflowVariable
        and:
        outs[0].inner.size() == 2
        and:
        outs[0].inner[0] instanceof CmdEvalParam
        outs[0].inner[0].getName() =~ /nxf_out_eval_\d+/
        (outs[0].inner[0] as CmdEvalParam).getTarget(binding) == 'this --one'
        and:
        outs[0].inner[1] instanceof CmdEvalParam
        outs[0].inner[1].getName() =~ /nxf_out_eval_\d+/
        (outs[0].inner[1] as CmdEvalParam).getTarget(binding) == 'tool --two'

    }
}
