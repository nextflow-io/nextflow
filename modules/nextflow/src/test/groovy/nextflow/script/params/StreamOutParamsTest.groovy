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
import nextflow.script.TestScriptRunner
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class StreamOutParamsTest extends Specification {

    def 'should create stream out param'() {

        setup:
        def text = '''
            process hola {
              output:
              stream z

              /echo ciao > $z/
            }
            '''

        def process = parseAndReturnProcess(text)

        when:
        def out0 = (StreamOutParam)process.config.getOutputs().get(0)

        then:
        out0.name == 'z'
        out0.outChannel instanceof DataflowQueue

    }

    def 'should replace stream topic' () {
        setup:
        def runner = new TestScriptRunner([process:[executor:'nope']])

        /*
         * Test a task with a very simple body.
         * For testing purposes the processor just return the script itself as result
         */
        when:
        def script =
                '''
            process sayHello  {
                output:
                stream z
                
                /echo Hello > $z/
            }
            '''

        runner.setScript(script).execute()

        // when no outputs are specified, the 'stdout' is the default output
        then:
        runner.result.val == "echo Hello world"
    }

}
