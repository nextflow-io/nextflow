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

package nextflow.script

import spock.lang.Specification

import java.nio.file.Paths

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import nextflow.processor.TaskProcessor
import static test.TestParser.parseAndReturnProcess
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ParamsInTest extends Specification {


    // ==============================================================
    //                  test *input* parameters
    // ==============================================================

    def testInputVals() {
        setup:
        def text = '''
            x = 'Hello'
            y = 'Hola'

            process hola {
              input:
              val x
              val x from y
              val x from 'ciao'
              val x from { [1,2] }

              return ''
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        def in1 = process.config.getInputs().get(0)
        def in2 = process.config.getInputs().get(1)
        def in3 = process.config.getInputs().get(2)
        def in4 = process.config.getInputs().get(3)

        then:
        process.config.getInputs().size() == 4

        in1.class == ValueInParam
        in1.name == 'x'
        in1.inChannel.val == 'Hello'

        in2.class == ValueInParam
        in2.name == 'x'
        in2.inChannel.val == 'Hola'

        in3.class == ValueInParam
        in3.name == 'x'
        in3.inChannel.val == 'ciao'

        in4.class == ValueInParam
        in4.name == 'x'
        in4.inChannel.val == 1
        in4.inChannel.val == 2
        in4.inChannel.val == Channel.STOP

    }

    def testFromMultiple() {


        setup:
        def text = '''
            A = 3
            B = 4

            process hola {
              input:
              val x from 1,2
              val y from ( {'a'}, {'b'} )
              val z from ( A, B )

              return ''
            }
            '''

        when:
        TaskProcessor process = parseAndReturnProcess(text)
        def in1 = process.config.getInputs().get(0)
        def in2 = process.config.getInputs().get(1)
        def in3 = process.config.getInputs().get(2)

        then:
        in1.name == 'x'
        in1.inChannel.val == 1
        in1.inChannel.val == 2
        in1.inChannel.val == Channel.STOP

        in2.name == 'y'
        in2.inChannel.val == 'a'
        in2.inChannel.val == 'b'
        in2.inChannel.val == Channel.STOP

        in3.name == 'z'
        in3.inChannel.val == 3
        in3.inChannel.val == 4
        in3.inChannel.val == Channel.STOP

    }

    def testFailFromMultipleChannels() {
        setup:
        def text = '''
            A = Channel.create()
            B = Channel.create()

            process hola {
              input:
              val x from A, B

              return ''
            }
            '''

        when:
        parseAndReturnProcess(text)
        then:
        thrown(IllegalArgumentException)

    }

    def testInputFiles() {
        setup:
        def text = '''
            x = java.nio.file.Paths.get('file.x')

            process hola {
              input:
              file x
              file f1 from x
              file f2 name 'abc' from x
              file f3:'*.fa' from x
              file 'file.txt' from x

              return ''
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        FileInParam in1 = process.config.getInputs().get(0)
        FileInParam in2 = process.config.getInputs().get(1)
        FileInParam in3 = process.config.getInputs().get(2)
        FileInParam in4 = process.config.getInputs().get(3)
        FileInParam in5 = process.config.getInputs().get(4)

        then:
        process.config.getInputs().size() == 5

        in1.name == 'x'
        in1.filePattern == '*'
        in1.inChannel.val == Paths.get('file.x')
        in1.index == 0

        in2.name == 'f1'
        in2.filePattern == '*'
        in2.inChannel.val == Paths.get('file.x')
        in2.index == 1

        in3.name == 'f2'
        in3.filePattern == 'abc'
        in3.inChannel.val == Paths.get('file.x')
        in3.index == 2

        in4.name == 'f3'
        in4.filePattern == '*.fa'
        in4.inChannel.val == Paths.get('file.x')
        in4.index == 3

        in5.name == 'file.txt'
        in5.filePattern == 'file.txt'
        in5.inChannel.val == Paths.get('file.x')
        in5.index == 4

    }

    def testInputFilesWithGString() {
        setup:
        def ctx = [x: 'main.txt', y: 'hello', z:'the_file_name']
        def text = '''
            q = java.nio.file.Paths.get('file.txt')

            process hola {
              input:
              file "$x" from q
              file "${y}.txt" from "str"
              file f2 name "${z}.fa" from q

              return ''
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        FileInParam in1 = process.config.getInputs().get(0)
        FileInParam in2 = process.config.getInputs().get(1)
        FileInParam in3 = process.config.getInputs().get(2)

        then:
        process.config.getInputs().size() == 3

        in1.name == '__$fileinparam<0>'
        in1.getFilePattern(ctx) == 'main.txt'
        in1.inChannel.val == Paths.get('file.txt')

        in2.name == '__$fileinparam<1>'
        in2.getFilePattern(ctx) == 'hello.txt'
        in2.inChannel.val == "str"

        in3.name == 'f2'
        in3.getFilePattern(ctx) == 'the_file_name.fa'
        in3.inChannel.val == Paths.get('file.txt')

    }

    def testInputFilesWithClosure() {
        setup:
        def ctx = [x: 'main.txt', y: 'hello', z:'the_file_name']
        def text = '''
            q = java.nio.file.Paths.get('file.txt')

            process hola {
              input:
              file { "$x" } from q
              file { "${y}.txt" } from "str"
              file f2 name { "${z}.fa" } from q
              file f3:{ "${z}.txt" } from q

              return ''
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        FileInParam in1 = process.config.getInputs().get(0)
        FileInParam in2 = process.config.getInputs().get(1)
        FileInParam in3 = process.config.getInputs().get(2)
        FileInParam in4 = process.config.getInputs().get(3)

        then:
        process.config.getInputs().size() == 4

        in1.name == '__$fileinparam<0>'
        in1.getFilePattern(ctx) == 'main.txt'
        in1.inChannel.val == Paths.get('file.txt')

        in2.name == '__$fileinparam<1>'
        in2.getFilePattern(ctx) == 'hello.txt'
        in2.inChannel.val == "str"

        in3.name == 'f2'
        in3.getFilePattern(ctx) == 'the_file_name.fa'
        in3.inChannel.val == Paths.get('file.txt')

        in4.name == 'f3'
        in4.getFilePattern(ctx) == 'the_file_name.txt'
        in4.inChannel.val == Paths.get('file.txt')

    }

    def testFromStdin() {
        setup:
        def text = '''
            x = 'Hola mundo'
            y = 'Ciao mondo'

            process hola {
              input:
              stdin x
              stdin from y

              return ''
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        def in1 = process.config.getInputs().get(0)
        def in2 = process.config.getInputs().get(1)

        then:
        process.config.getInputs().size() == 2

        in1.class == StdInParam
        in1.name == '-'
        in1.inChannel.val == 'Hola mundo'

        in2.class == StdInParam
        in2.name == '-'
        in2.inChannel.val == 'Ciao mondo'
    }

    def testInputEnv() {
        setup:
        def text = '''
            x = 'aaa'
            y = [1,2]

            process hola {
              input:
              env VAR_X from x
              env 'VAR_Y' from y

              return ''
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        def in1 = process.config.getInputs().get(0)
        def in2 = process.config.getInputs().get(1)

        then:
        process.config.getInputs().size() == 2

        in1.class == EnvInParam
        in1.name == 'VAR_X'
        in1.inChannel.val == 'aaa'

        in2.class == EnvInParam
        in2.name == 'VAR_Y'
        in2.inChannel.val == 1
        in2.inChannel.val == 2
        in2.inChannel.val == Channel.STOP

    }


    def testInputMap() {
        setup:
        def text = '''
            x = 'Hola mundo'

            process hola {
              input:
              set (p) from x
              set (p, q) from x
              set (v, 'file_name.fa') from 'str'
              set (p, 'file_name.txt', '-' ) from { 'ciao' }
              set (t, [file: 'file.fa'] ) from 0

              return ''
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        SetInParam in1 = process.config.getInputs().get(0)
        SetInParam in2 = process.config.getInputs().get(1)
        SetInParam in3 = process.config.getInputs().get(2)
        SetInParam in4 = process.config.getInputs().get(3)
        SetInParam in5 = process.config.getInputs().get(4)

        then:
        process.config.getInputs().size() == 5

        in1.inner.size() == 1
        in1.inner.get(0) instanceof ValueInParam
        in1.inner.get(0).index == 0
        in1.inner.get(0).mapIndex == 0
        in1.inner.get(0).name == 'p'
        in1.inChannel.val == 'Hola mundo'

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

        in5.inner.size() == 2
        in5.inner.get(0) instanceof ValueInParam
        in5.inner.get(0).name == 't'
        in5.inner.get(0).index == 4
        in5.inner.get(0).mapIndex == 0
        in5.inner.get(1) instanceof FileInParam
        in5.inner.get(1).name == 'file'
        in5.inner.get(1).filePattern == 'file.fa'
        in5.inner.get(1).index == 4
        in5.inner.get(1).mapIndex == 1
        in5.inChannel.val == 0

    }

    def testSetFileWithGString() {

        setup:
        def text = '''
            q = 'the file content'

            process hola {
              input:
              set( 'name_$x' ) from q
              set( "${x}_name.${str}" ) from q

              set( file("hola_${x}") ) from q
              set( file( handle: "${x}.txt") ) from q

              set file( { "${x}_name.txt" } ) from q
              set file( handle: { "name_${x}.txt" } ) from q

              return ''
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        SetInParam in0 = process.config.getInputs().get(0)
        SetInParam in1 = process.config.getInputs().get(1)
        SetInParam in2 = process.config.getInputs().get(2)
        SetInParam in3 = process.config.getInputs().get(3)
        SetInParam in4 = process.config.getInputs().get(4)
        SetInParam in5 = process.config.getInputs().get(5)
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
        (in3.inner[0] as FileInParam).name == 'handle'
        (in3.inner[0] as FileInParam).getFilePattern(ctx) == 'the_file.txt'

        in4.inner[0] instanceof FileInParam
        (in4.inner[0] as FileInParam).name == '__$fileinparam<4:0>'
        (in4.inner[0] as FileInParam).getFilePattern(ctx) == 'the_file_name.txt'

        in5.inner[0] instanceof FileInParam
        (in5.inner[0] as FileInParam).name == 'handle'
        (in5.inner[0] as FileInParam).getFilePattern(ctx) == 'name_the_file.txt'
    }

    def testInputMap2() {
        setup:
        def text = '''
            process hola {
              input:
              set( val(a), file(x), val(b) ) from 1
              set( val(p), file('txt'), env('q') ) from 2
              set( val(v), file(xx:'yy'), stdin, env(W) ) from 3

              return ''
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        SetInParam in0 = process.config.getInputs().get(0)
        SetInParam in1 = process.config.getInputs().get(1)
        SetInParam in2 = process.config.getInputs().get(2)

        then:
        process.config.getInputs().size() == 3

        in0.name == '__$setinparam<0>'
        in0.inChannel.val == 1
        in0.inner.size() == 3
        in0.inner.get(0) instanceof ValueInParam
        in0.inner.get(0).name == 'a'
        in0.inner.get(0).mapIndex == 0
        in0.inner.get(1) instanceof FileInParam
        in0.inner.get(1).name == 'x'
        in0.inner.get(1).filePattern == '*'
        in0.inner.get(1).mapIndex == 1
        in0.inner.get(2) instanceof ValueInParam
        in0.inner.get(2).name == 'b'
        in0.inner.get(2).mapIndex == 2

        in1.inner.get(0) instanceof ValueInParam
        in1.inner.get(0).name == 'p'
        in1.inner.get(0).mapIndex == 0
        in1.inner.get(1) instanceof FileInParam
        in1.inner.get(1).name == 'txt'
        in1.inner.get(1).filePattern == 'txt'
        in1.inner.get(1).mapIndex == 1
        in1.inner.get(2) instanceof EnvInParam
        in1.inner.get(2).name == 'q'
        in1.inner.get(2).mapIndex == 2

        in2.inner.get(0) instanceof ValueInParam
        in2.inner.get(0).name == 'v'
        in2.inner.get(0).mapIndex == 0
        in2.inner.get(1) instanceof FileInParam
        in2.inner.get(1).name == 'xx'
        in2.inner.get(1).filePattern == 'yy'
        in2.inner.get(1).mapIndex == 1
        in2.inner.get(2) instanceof StdInParam
        in2.inner.get(2).name == '-'
        in2.inner.get(2).mapIndex == 2
        in2.inner.get(3) instanceof EnvInParam
        in2.inner.get(3).name == 'W'
        in2.inner.get(3).mapIndex == 3

    }


    def testEachInParam() {

        setup:
        def text = '''
            x = 'aaa'
            y = [1,2]

            process hola {
              input:
              each x
              each p from y
              each z from q
              each file(foo) from foo_ch
              each file('bar') from bar_ch

              return ''
            }
            '''
        when:

        def binding =  [:]
        binding.q = new DataflowQueue<>()
        binding.q << 1 << 2 << 3  << Channel.STOP
        binding.foo_ch = 'file-a.txt'
        binding.bar_ch = 'file-x.fa'

        def process = parseAndReturnProcess(text, binding)
        def in0 = (EachInParam)process.config.getInputs().get(0)
        def in1 = (EachInParam)process.config.getInputs().get(1)
        def in2 = (EachInParam)process.config.getInputs().get(2)
        def in3 = (EachInParam)process.config.getInputs().get(3)
        def in4 = (EachInParam)process.config.getInputs().get(4)

        then:
        process.config.getInputs().size() == 5

        in0.class == EachInParam
        in0.inChannel instanceof DataflowVariable
        in0.inChannel.val == ['aaa']
        in0.inner.name == 'x'
        in0.inner.owner == in0

        in1.class == EachInParam
        in1.name == '__$eachinparam<1>'
        in1.inChannel instanceof DataflowVariable
        in1.inChannel.val == [1,2]
        in1.inner.name == 'p'
        in1.inner instanceof ValueInParam
        in1.inner.owner == in1

        in2.class == EachInParam
        in2.name == '__$eachinparam<2>'
        in2.inChannel.val == [1,2,3]
        in2.inner instanceof ValueInParam
        in2.inner.name == 'z'
        in2.inner.owner == in2

        in3.class == EachInParam
        in3.name == '__$eachinparam<3>'
        in3.inChannel instanceof DataflowVariable
        in3.inChannel.val == ['file-a.txt']
        in3.inner instanceof FileInParam
        in3.inner.name == 'foo'
        in3.inner.owner == in3

        in4.class == EachInParam
        in4.name == '__$eachinparam<4>'
        in4.inChannel instanceof DataflowVariable
        in4.inChannel.val == ['file-x.fa']
        in4.inner instanceof FileInParam
        in4.inner.name == 'bar'
        in4.inner.filePattern == 'bar'
        in4.inner.owner == in4

    }

    def testMissingVariable () {

        setup:
        def text = '''

            process hola {
              input:
                val x

              return ''
            }
            '''
        when:
        def process = parseAndReturnProcess(text)
        process.config.getInputs().get(0).inChannel

        then:
        thrown(MissingPropertyException)

    }


    def 'should decode param inputs ' () {

        def param
        def holder = []

        when:
        param = new ValueInParam(Mock(Binding), holder)
        then:
        param.decodeInputs( ['a','b','c'] ) == 'a'

        when:
        param = new ValueInParam(Mock(Binding), holder)
        then:
        param.decodeInputs( ['a','b','c'] ) == 'b'

        when:
        param = new ValueInParam(Mock(Binding), [])
        param.owner = new EachInParam(Mock(Binding), [])
        then:
        param.decodeInputs( ['a','b','c'] ) == 'a'

        when:
        param = new ValueInParam(Mock(Binding), [])
        param.owner = new EachInParam(Mock(Binding), [])
        then:
        param.decodeInputs( [[1,2,3],'b','c'] ) == [1,2,3]
    }

    def testInvalidInputs() {

        given:
        TaskProcessor process
        Exception e

        when:
        process = parseAndReturnProcess('''
                        x = null
                        process hola {
                          input:
                            val x
            
                          return ''
                        }
                        ''')
        process.config.getInputs().get(0).inChannel
        then:
        e = thrown(IllegalArgumentException)
        e.message == 'A process input channel evaluates to null -- Invalid declaration `val x`'

        when:
        process = parseAndReturnProcess('''
                        x = null
                        process hola {
                          input:
                            file foo from x
            
                          return ''
                        }
                        ''')
        process.config.getInputs().get(0).inChannel
        then:
        e = thrown(IllegalArgumentException)
        e.message == 'A process input channel evaluates to null -- Invalid declaration `file foo`'

    }

    def testInterpolateGlobalVars() {

        given:
        def TEXT = '''
            foo = 1
            bar = 2
            baz = 3 
            alpha = [beta: 'hello']
            
            process foo {
              input:
              val x from "$foo"
              val y from "${alpha.beta}"
              val x from Channel.value("$bar")
              val z from Channel.from("$baz")
              
              script:
              """
              echo $x
              """
            }
            '''
        
        when:
        def process = parseAndReturnProcess(TEXT)
        def in0 = (ValueInParam)process.config.getInputs().get(0)
        def in1 = (ValueInParam)process.config.getInputs().get(1)
        def in2 = (ValueInParam)process.config.getInputs().get(2)
        def in3 = (ValueInParam)process.config.getInputs().get(3)
        then:
        in0.inChannel.val == "1"
        in1.inChannel.val == "hello"
        in2.inChannel.val == '2'
        in3.inChannel.val == '3'


    }

}
