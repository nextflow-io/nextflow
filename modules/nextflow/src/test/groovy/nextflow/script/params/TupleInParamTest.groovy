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
class TupleInParamTest extends Specification {

    def 'should create input tuples'() {
        setup:
        def text = '''
            x = 'Hola mundo'

            process hola {
              input:
              tuple val(p) from x
              tuple val(p), val(q) from x
              tuple val(v), path('file_name.fa') from 'str'
              tuple val(p), path('file_name.txt'), '-' from { 'ciao' }
              tuple val(p), path(z, stageAs: 'file*')
              
              /foo/
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        TupleInParam in1 = process.config.getInputs().get(0)
        TupleInParam in2 = process.config.getInputs().get(1)
        TupleInParam in3 = process.config.getInputs().get(2)
        TupleInParam in4 = process.config.getInputs().get(3)
        TupleInParam in5 = process.config.getInputs().get(4)

        then:
        in1.inner.size() == 1
        in1.inner.get(0) instanceof ValueInParam
        in1.inner.get(0).index == 0
        in1.inner.get(0).mapIndex == 0
        in1.inner.get(0).name == 'p'
        in1.inChannel.val == 'Hola mundo'
        and:
        in2.inner.size() == 2
        in2.inner.get(0) instanceof ValueInParam
        in2.inner.get(0).name == 'p'
        in2.inner.get(0).index == 1
        in2.inner.get(0).mapIndex == 0
        in2.inner.get(1) instanceof ValueInParam
        in2.inner.get(1).name == 'q'
        in2.inner.get(1).index == 1
        in2.inner.get(1).mapIndex == 1
        in2.inChannel.val == 'Hola mundo'
        and:
        in3.inner.size() == 2
        in3.inner.get(0) instanceof ValueInParam
        in3.inner.get(0).name == 'v'
        in3.inner.get(0).index == 2
        in3.inner.get(0).mapIndex == 0
        in3.inner.get(1) instanceof FileInParam
        in3.inner.get(1).name == 'file_name.fa'
        in3.inner.get(1).filePattern == 'file_name.fa'
        in3.inner.get(1).index == 2
        in3.inner.get(1).mapIndex == 1
        in3.inChannel.val == 'str'
        and:
        in4.inner.size() == 3
        in4.inner.get(0) instanceof ValueInParam
        in4.inner.get(0).name == 'p'
        in4.inner.get(0).index == 3
        in4.inner.get(0).mapIndex == 0
        in4.inner.get(1) instanceof FileInParam
        in4.inner.get(1).name == 'file_name.txt'
        in4.inner.get(1).filePattern == 'file_name.txt'
        in4.inner.get(1).index == 3
        in4.inner.get(1).mapIndex == 1
        in4.inner.get(2) instanceof StdInParam
        in4.inner.get(2).name == '-'
        in4.inner.get(2).index == 3
        in4.inner.get(2).mapIndex == 2
        in4.inChannel.val == 'ciao'
        and:
        in5.inner.size() == 2
        in5.inner.get(0) instanceof ValueInParam
        in5.inner.get(0).index == 4
        in5.inner.get(0).mapIndex == 0
        in5.inner.get(1) instanceof FileInParam
        in5.inner.get(1).name == 'z'
        in5.inner.get(1).filePattern == 'file*'
        in5.inner.get(1).index == 4
        in5.inner.get(1).mapIndex == 1

    }

    def 'should create tuple with gstring'() {

        setup:
        def text = '''
            q = 'the file content'

            process hola {
              input:
              tuple( 'name_$x' ) from q
              tuple( "${x}_name.${str}" ) from q

              tuple( file("hola_${x}") ) from q

              tuple file( { "${x}_name.txt" } ) from q

              /foo/
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        TupleInParam in0 = process.config.getInputs().get(0)
        TupleInParam in1 = process.config.getInputs().get(1)
        TupleInParam in2 = process.config.getInputs().get(2)
        TupleInParam in3 = process.config.getInputs().get(3)
        def ctx = [x:'the_file', str: 'fastq']

        then:
        in0.inChannel.val == 'the file content'
        in0.inner[0] instanceof FileInParam
        (in0.inner[0] as FileInParam).name == 'name_$x'
        (in0.inner[0] as FileInParam).getFilePattern(ctx) == 'name_$x'

        in1.inner[0] instanceof FileInParam
        (in1.inner[0] as FileInParam).name == '__$fileinparam<1:0>'
        (in1.inner[0] as FileInParam).getFilePattern(ctx) == 'the_file_name.fastq'

        in2.inner[0] instanceof FileInParam
        (in2.inner[0] as FileInParam).name == '__$fileinparam<2:0>'
        (in2.inner[0] as FileInParam).getFilePattern(ctx) == 'hola_the_file'

        in3.inner[0] instanceof FileInParam
        (in3.inner[0] as FileInParam).name == '__$fileinparam<3:0>'
        (in3.inner[0] as FileInParam).getFilePattern(ctx) == 'the_file_name.txt'
    }


}
