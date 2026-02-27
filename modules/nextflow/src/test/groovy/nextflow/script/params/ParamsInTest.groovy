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

import java.nio.file.Paths

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.*
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class ParamsInTest extends Dsl2Spec {

    // ==============================================================
    //                  test *input* parameters
    // ==============================================================

    def testInputVals() {
        setup:
        def text = '''
            process hola {
              input:
              val w
              val x
              val y
              val z

              script:
              ''
            }

            workflow {
              hola('Hello', 'Hola', 'ciao', channel.of(1, 2))
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
        in1.name == 'w'
        in1.inChannel.val == 'Hello'

        in2.class == ValueInParam
        in2.name == 'x'
        in2.inChannel.val == 'Hola'

        in3.class == ValueInParam
        in3.name == 'y'
        in3.inChannel.val == 'ciao'

        in4.class == ValueInParam
        in4.name == 'z'
        in4.inChannel.val == 1
        in4.inChannel.val == 2
        in4.inChannel.val == Channel.STOP

    }

    def testFromMultiple() {

        setup:
        def text = '''
            process hola {
              input:
              val x
              val y
              val z

              script:
              ''
            }

            workflow {
              A = 3
              B = 4

              def x = channel.of(1,2)
              def y = channel.of('a', 'b')
              def z = channel.of(A, B)
              hola(x, y, z)
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
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

    def testInputFiles() {
        setup:
        def text = '''
            process hola {
              input:
              file x
              file f1
              file 'file.txt'

              script:
              ''
            }

            workflow {
              x = file('file.x')

              hola(x, x, x)
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        FileInParam in1 = process.config.getInputs().get(0)
        FileInParam in2 = process.config.getInputs().get(1)
        FileInParam in3 = process.config.getInputs().get(2)

        then:
        process.config.getInputs().size() == 3

        in1.name == 'x'
        in1.filePattern == '*'
        in1.inChannel.val.name == 'file.x'
        in1.index == 0

        in2.name == 'f1'
        in2.filePattern == '*'
        in2.inChannel.val.name == 'file.x'
        in2.index == 1

        in3.name == 'file.txt'
        in3.filePattern == 'file.txt'
        in3.inChannel.val.name == 'file.x'
        in3.index == 2

    }

    def testInputFilesWithGString() {
        setup:
        def ctx = [id: 'hello']
        def text = '''
            process hola {
              input:
              val id
              file "${id}.txt"

              script:
              ''
            }

            workflow {
              hola('hello', 'str')
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        def in1 = process.config.getInputs().get(0)
        def in2 = process.config.getInputs().get(1)

        then:
        process.config.getInputs().size() == 2

        in1.name == 'id'
        in1.inChannel.val == 'hello'

        in2.name == '__$fileinparam<1>'
        in2.getFilePattern(ctx) == 'hello.txt'
        in2.inChannel.val == "str"

    }

    def testFromStdin() {
        setup:
        def text = '''
            process hola {
              input:
              stdin
              stdin

              script:
              ''
            }

            workflow {
              x = 'Hola mundo'
              y = 'Ciao mondo'

              hola(x, y)
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
            process hola {
              input:
              env 'VAR_X'
              env 'VAR_Y'

              script:
              ''
            }

            workflow {
              x = 'aaa'
              y = channel.of(1,2)

              hola(x, y)
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

    def testInputTuple2() {
        setup:
        def text = '''
            process hola {
              input:
              tuple( val(a), file(x), val(b) )
              tuple( val(p), file('txt'), env('q') )
              tuple( val(v), file(xx:'yy'), stdin, env('W') )

              script:
              ''
            }

            workflow {
              hola(1, 2, 3)
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        TupleInParam in0 = process.config.getInputs().get(0)
        TupleInParam in1 = process.config.getInputs().get(1)
        TupleInParam in2 = process.config.getInputs().get(2)

        then:
        process.config.getInputs().size() == 3

        in0.name == '__$tupleinparam<0>'
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
            process hola {
              input:
              each x
              each p
              each z
              each file('foo')
              each file('bar')

              script:
              ''
            }

            workflow {
              hola('aaa', [1, 2], channel.of(1, 2, 3), 'file-a.txt', 'file-x.fa')
            }
            '''
        when:
        def process = parseAndReturnProcess(text)
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

    /*
     * test path qualifier
     */

    def 'test input paths'() {
        setup:
        def FILE = '/data/file.txt'
        def text = """
            process hola {
              input:
              path x, arity: '1'
              path f1, arity: '1..2'
              path '*.fa', arity: '1..*'
              path 'file.txt'
              path f2, name: '*.fa'
              path f3, stageAs: '*.txt'

              script:
              ''
            }

            workflow {
              x = '$FILE'
              hola(x, x, x, x, x, x)
            }
            """

        when:
        def process = parseAndReturnProcess(text)
        FileInParam in0 = process.config.getInputs().get(0)
        FileInParam in1 = process.config.getInputs().get(1)
        FileInParam in2 = process.config.getInputs().get(2)
        FileInParam in3 = process.config.getInputs().get(3)
        FileInParam in4 = process.config.getInputs().get(4)
        FileInParam in5 = process.config.getInputs().get(5)

        then:

        in0.name == 'x'
        in0.filePattern == '*'
        in0.inChannel.val == FILE
        in0.index == 0
        in0.isPathQualifier()
        in0.arity == new ArityParam.Range(1, 1)

        in1.name == 'f1'
        in1.filePattern == '*'
        in1.inChannel.val == FILE
        in1.index == 1
        in1.isPathQualifier()
        in1.arity == new ArityParam.Range(1, 2)

        in2.name == '*.fa'
        in2.filePattern == '*.fa'
        in2.inChannel.val == FILE
        in2.index == 2
        in2.isPathQualifier()
        in2.arity == new ArityParam.Range(1, Integer.MAX_VALUE)

        in3.name == 'file.txt'
        in3.filePattern == 'file.txt'
        in3.inChannel.val == FILE
        in3.index == 3
        in3.isPathQualifier()

        in4.name == 'f2'
        in4.filePattern == '*.fa'
        in4.index == 4
        in4.isPathQualifier()

        in5.name == 'f3'
        in5.filePattern == '*.txt'
        in5.index == 5
        in5.isPathQualifier()
    }

    def 'test input paths with gstring'() {
        setup:
        def ctx = [id: 'hello']
        def text = '''
            process hola {
              input:
              val id
              path "${id}.txt"

              script:
              ''
            }

            workflow {
              hola('hello', 'str')
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        def in1 = process.config.getInputs().get(0)
        def in2 = process.config.getInputs().get(1)

        then:
        process.config.getInputs().size() == 2

        in1.name == 'id'
        in1.inChannel.val == 'hello'

        in2.name == '__$pathinparam<1>'
        in2.getFilePattern(ctx) == 'hello.txt'
        in2.inChannel.val == "str"
        in2.isPathQualifier()
    }

    def 'test tuple path with gstring'() {

        setup:
        def text = '''
            process hola {
              input:
              tuple val(id), path("hola_${id}")

              script:
              ''
            }

            workflow {
              hola(['the_file', '/the/file/path'])
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        TupleInParam in1 = process.config.getInputs().get(0)
        def ctx = [id:'the_file']

        then:
        in1.inChannel.val == ['the_file', '/the/file/path']
        in1.inner[0] instanceof ValueInParam
        (in1.inner[0] as ValueInParam).getName() == 'id'

        in1.inner[1] instanceof FileInParam
        (in1.inner[1] as FileInParam).name == '__$pathinparam<0:1>'
        (in1.inner[1] as FileInParam).getFilePattern(ctx) == 'hola_the_file'
        (in1.inner[1] as FileInParam).isPathQualifier()

    }

    def 'test input path with map'() {
        setup:
        def text = '''
            process hola {
              input:
              tuple( val(a), path(x) )
              tuple( val(p), path('txt') )
              tuple( val(v), path(xx, stageAs: 'yy') )

              script:
              ''
            }

            workflow {
                hola(1,2,3)
            }
            '''

        when:
        def process = parseAndReturnProcess(text)
        TupleInParam in0 = process.config.getInputs().get(0)
        TupleInParam in1 = process.config.getInputs().get(1)
        TupleInParam in2 = process.config.getInputs().get(2)

        then:
        process.config.getInputs().size() == 3

        in0.name == '__$tupleinparam<0>'
        in0.inChannel.val == 1
        in0.inner.size() == 2
        in0.inner.get(0) instanceof ValueInParam
        in0.inner.get(0).name == 'a'
        in0.inner.get(0).mapIndex == 0
        in0.inner.get(1) instanceof FileInParam
        in0.inner.get(1).name == 'x'
        in0.inner.get(1).filePattern == '*'
        in0.inner.get(1).mapIndex == 1
        in0.inner.get(1).isPathQualifier()

        in1.inner.get(0) instanceof ValueInParam
        in1.inner.get(0).name == 'p'
        in1.inner.get(0).mapIndex == 0
        in1.inner.get(1) instanceof FileInParam
        in1.inner.get(1).name == 'txt'
        in1.inner.get(1).filePattern == 'txt'
        in1.inner.get(1).mapIndex == 1
        in1.inner.get(1).isPathQualifier()

        in2.inner.get(0) instanceof ValueInParam
        in2.inner.get(0).name == 'v'
        in2.inner.get(0).mapIndex == 0
        in2.inner.get(1) instanceof FileInParam
        in2.inner.get(1).name == 'xx'
        in2.inner.get(1).filePattern == 'yy'
        in2.inner.get(1).mapIndex == 1
        in2.inner.get(1).isPathQualifier()

    }


    def 'test input each path'() {

        setup:
        def text = '''

            process hola {
              input:
              each path('foo')
              each path('bar')

              script:
              ''
            }

            workflow {
              hola('file-a.txt', 'file-x.fa')
            }
            '''
        when:

        def binding =  [:]

        def process = parseAndReturnProcess(text, binding)
        def in0 = (EachInParam)process.config.getInputs().get(0)
        def in1 = (EachInParam)process.config.getInputs().get(1)

        then:
        process.config.getInputs().size() == 2

        in0.class == EachInParam
        in0.name == '__$eachinparam<0>'
        in0.inChannel instanceof DataflowVariable
        in0.inChannel.val == ['file-a.txt']
        in0.inner instanceof FileInParam
        (in0.inner as FileInParam).name == 'foo'
        (in0.inner as FileInParam).owner == in0
        (in0.inner as FileInParam).isPathQualifier()

        in1.class == EachInParam
        in1.name == '__$eachinparam<1>'
        in1.inChannel instanceof DataflowVariable
        in1.inChannel.val == ['file-x.fa']
        in1.inner instanceof FileInParam
        (in1.inner as FileInParam).name == 'bar'
        (in1.inner as FileInParam).filePattern == 'bar'
        (in1.inner as FileInParam).owner == in1
        (in1.inner as FileInParam).isPathQualifier()

    }


    def 'should check is tuple item' () {

        setup:
        def text = '''
            process hola {
              input:
              val x
              tuple val(y), file(z)

              script:
              /command/
            }

            workflow {
              ch = 'something'
              hola(ch, ch)
            }
            '''
        when:


        def process = parseAndReturnProcess(text)
        def in0 = (ValueInParam)process.config.getInputs().get(0)
        def in1 = (TupleInParam)process.config.getInputs().get(1)

        then:
        !in0.isNestedParam()
        !in1.isNestedParam()
        (in1.inner[0] as ValueInParam).isNestedParam()
        (in1.inner[1] as FileInParam).isNestedParam()
    }

}
