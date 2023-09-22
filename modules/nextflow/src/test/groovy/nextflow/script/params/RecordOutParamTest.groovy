/*
 * Copyright 2023, Seqera Labs
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

import groovyx.gpars.dataflow.DataflowVariable
import test.Dsl2Spec
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class RecordOutParamTest extends Dsl2Spec {

    def 'should define output maps'() {

        setup:
        def text = '''
            process hola {
              output:
                record( x: val(x) )
                record( y: val(y), stdout: stdout, fa: file('*.fa')  )
                record( stdout: stdout, z: val(z) )
                record( foo: env(FOO), bar: env(BAR) )
              return ''
            }
            
            workflow {
              hola()
            }
            '''

        def process = parseAndReturnProcess(text)

        when:
        RecordOutParam out0 = process.config.getOutputs().get(0)
        RecordOutParam out1 = process.config.getOutputs().get(1)
        RecordOutParam out2 = process.config.getOutputs().get(2)
        RecordOutParam out3 = process.config.getOutputs().get(3)

        then:
        process.config.getOutputs().size() == 4

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
        out1.inner.size() == 3

        out2.outChannel instanceof DataflowVariable
        out2.inner.size() == 2
        out2.inner[0] instanceof StdOutParam
        out2.inner[0].name == '-'
        out2.inner[0].index == 2
        out2.inner[1] instanceof ValueOutParam
        out2.inner[1].name == 'z'
        out2.inner[1].index == 2

        out3.outChannel instanceof DataflowVariable
        out3.inner.size() ==2
        out3.inner[0] instanceof EnvOutParam
        out3.inner[0].getName() == 'FOO'
        out3.inner[1] instanceof EnvOutParam
        out3.inner[1].getName() == 'BAR'

    }

}