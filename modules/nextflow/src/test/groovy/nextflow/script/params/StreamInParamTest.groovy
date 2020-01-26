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

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class StreamInParamTest extends Specification {

    def 'should create stream in param'() {
        setup:
        def text = '''
            process foo {
              input: 
              stream x  
                
              /cat < $x/
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        StreamInParam in1 = process.config.getInputs().get(0)

        then:
        process.config.getInputs().size() == 1

        in1.name == 'x'
//        in1.inChannel.val == Paths.get('file.x')
//        in1.index == 0


    }
}
