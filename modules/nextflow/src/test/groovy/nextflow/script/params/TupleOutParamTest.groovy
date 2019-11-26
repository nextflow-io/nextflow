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

import static test.TestParser.*

import groovyx.gpars.dataflow.DataflowQueue
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TupleOutParamTest extends Specification {

    def 'should define output tuples'() {

        setup:
        def text = '''
            process hola {
              output:
                tuple(x) into p
                tuple(y, '-', '*.fa') into q mode flatten
                tuple(stdout, z) into t mode combine

              return ''
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

        out0.outChannel instanceof DataflowQueue
        out0.outChannel == binding.p
        out0.inner.size() == 1
        out0.inner[0] instanceof ValueOutParam
        out0.inner[0].name == 'x'
        out0.inner[0].index == 0
        out0.mode == BasicMode.standard

        out1.outChannel instanceof DataflowQueue
        out1.outChannel == binding.q
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
        out1.mode == BasicMode.flatten

        out2.outChannel instanceof DataflowQueue
        out2.outChannel == binding.t
        out2.inner.size() == 2
        out2.inner[0] instanceof StdOutParam
        out2.inner[0].name == '-'
        out2.inner[0].index == 2
        out2.inner[1] instanceof ValueOutParam
        out2.inner[1].name == 'z'
        out2.inner[1].index == 2
        out2.mode == TupleOutParam.CombineMode.combine

    }

    def 'should define output tuples 2'() {

        setup:
        def text = '''
            process hola {
              output:
                tuple val(x) into p
                tuple val(y), stdout, file('*.fa') into q mode flatten
                tuple stdout, val(z) into t mode combine

              return ''
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

        out0.outChannel instanceof DataflowQueue
        out0.outChannel == binding.p
        out0.inner.size() == 1
        out0.inner[0] instanceof ValueOutParam
        out0.inner[0].name == 'x'
        out0.inner[0].index == 0
        out0.mode == BasicMode.standard

        out1.outChannel instanceof DataflowQueue
        out1.outChannel == binding.q
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
        out1.mode == BasicMode.flatten

        out2.outChannel instanceof DataflowQueue
        out2.outChannel == binding.t
        out2.inner.size() == 2
        out2.inner[0] instanceof StdOutParam
        out2.inner[0].name == '-'
        out2.inner[0].index == 2
        out2.inner[1] instanceof ValueOutParam
        out2.inner[1].name == 'z'
        out2.inner[1].index == 2
        out2.mode == TupleOutParam.CombineMode.combine

    }

    def 'should create tuple of env' () {
        setup:
        def text = '''
            process hola {
              output:
                tuple env(FOO), env(BAR) into ch
              
              /echo command/ 
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
        outs[0].outChannel instanceof DataflowQueue
        outs[0].outChannel == binding.ch
        and:
        outs[0].inner.size() ==2
        and:
        outs[0].inner[0] instanceof EnvOutParam
        outs[0].inner[0].getName() == 'FOO'
        and:
        outs[0].inner[1] instanceof EnvOutParam
        outs[0].inner[1].getName() == 'BAR'
    }
}
