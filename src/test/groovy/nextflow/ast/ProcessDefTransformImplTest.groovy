/*
 * Copyright (c) 2012, the authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.ast

import java.nio.file.Paths

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import nextflow.Nextflow
import nextflow.script.EachInParam
import nextflow.script.EnvInParam
import nextflow.script.FileInParam
import nextflow.script.FileOutParam
import nextflow.script.FileSharedParam
import nextflow.script.StdInParam
import nextflow.script.StdOutParam
import nextflow.processor.TaskProcessor
import nextflow.script.ValueInParam
import nextflow.script.ValueOutParam
import nextflow.script.ValueSharedParam
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessDefTransformImplTest extends Specification {

    def parse ( String scriptText, Map map = [:] ) {
        new TestParser().parse(scriptText, map)
    }


    // ==============================================================
    //                  test *input* parameters
    // ==============================================================

    def testInputVal() {

        setup:
        def text = '''
            process hola {
              input:
              val x

              return ''
            }
            '''

        // when 'x' variable is not declared in the script context
        // throws an exception
        when:
        parse(text).run()
        then:
        thrown(MissingPropertyException)

        when:
        TaskProcessor process = parse(text,[x: 1]).run()
        ValueInParam in1 = process.taskConfig.getInputs().get(0)
        then:
        in1.name == 'x'
        in1.inTarget == 1
        in1.getInChannel() instanceof DataflowVariable
        in1.getInChannel().getVal() == 1

    }

    def testInputVal2() {

        setup:
        def text = '''
            process hola {
              input:
              val 'hola' as q
              val x as p1
              val y as p2

              return ''
            }
            '''

        def y = new DataflowQueue()

        when:
        TaskProcessor process = parse(text,[x: 1, y:y ]).run()
        ValueInParam in1 = process.taskConfig.getInputs().get(0)
        ValueInParam in2 = process.taskConfig.getInputs().get(1)
        ValueInParam in3 = process.taskConfig.getInputs().get(2)

        then:
        in1.name == 'q'
        in1.inTarget == 'hola'
        in1.inChannel instanceof DataflowVariable
        in1.inChannel.getVal() == 'hola'

        in2.name == 'p1'
        in2.inTarget == 1
        in2.inChannel instanceof DataflowVariable
        in2.inChannel.getVal() == 1

        in3.name == 'p2'
        in3.inTarget == y
        in3.inChannel instanceof DataflowQueue
        in3.inChannel  == y

    }

    def testInputVal3() {

        setup:
        def text = '''
            X = 1
            process hola {
              input:
              val X
              val 2 as y
              val { 'str' } as w
              val ([1,2,3]) as z

              return ''
            }
            '''

        when:
        TaskProcessor process = parse(text).run()
        ValueInParam in1 = process.taskConfig.getInputs().get(0)
        ValueInParam in2 = process.taskConfig.getInputs().get(1)
        ValueInParam in3 = process.taskConfig.getInputs().get(2)
        ValueInParam in4 = process.taskConfig.getInputs().get(3)

        then:
        in1.name == 'X'
        in1.inTarget == 1
        in1.inChannel instanceof DataflowVariable
        in1.inChannel.val == 1

        in2.name == 'y'
        in2.inTarget == 2
        in2.inChannel instanceof DataflowVariable
        in2.inChannel.val == 2

        in3.name == 'w'
        in3.inTarget instanceof Closure
        in3.inChannel instanceof DataflowVariable
        in3.inChannel.val == 'str'

        in4.name == 'z'
        in4.inTarget == [1,2,3]
        in4.inChannel instanceof DataflowQueue
        in4.inChannel.val == 1
        in4.inChannel.val == 2
        in4.inChannel.val == 3
        in4.inChannel.val == Channel.STOP

    }


    def testInputFile() {

        setup:
        def text = '''
            blast_x = 'blah blah'
            blast_y = 'blah blah'
            blast_z = 'blah blah'

            process hola {
              input:
              file blast_x
              file blast_y as file_xxx
              file blast_z as 'file.txt'
              file { 1 }
              file 'hello' as var_name as 'file.name'

              return ''
            }
            '''


        when:
        TaskProcessor process = parse(text).run()
        FileInParam in1 = process.taskConfig.getInputs().get(0)
        FileInParam in2 = process.taskConfig.getInputs().get(1)
        FileInParam in3 = process.taskConfig.getInputs().get(2)
        FileInParam in4 = process.taskConfig.getInputs().get(3)
        FileInParam in5 = process.taskConfig.getInputs().get(4)

        then:
        // verifies:
        // - the variable name is 'blast_x'
        // - the filePattern is null
        // - the channel returns the variable value as defined in the script
        in1.name == 'blast_x'
        in1.filePattern == null
        in1.inChannel.getVal() == 'blah blah'

        // - the name is the one defined by the 'as' keyword
        // - the filePattern is undefined
        // = the channel binds to 'blast_y' value
        in2.name == 'file_xxx'
        in2.filePattern == null
        in2.inChannel.getVal() == 'blah blah'

        // - the (file) name get the value specified by the 'as' keyword
        // - target get the 'blast_z' value
        // = the channel binds to 'blast_z' value
        in3.name == 'blast_z'
        in3.filePattern == 'file.txt'
        in3.inChannel.getVal() == 'blah blah'

        // the file argument is not a variable, but an expression, thus
        // - the (file) name is just '*'
        // - the target is the closure
        // - the channel holds the closure value
        in4.name == null
        in4.filePattern == null
        in4.inChannel.getVal() == 1

        in5.name == 'var_name'
        in5.filePattern == 'file.name'
        in5.inChannel.getVal() == 'hello'

    }

    def tesInputEnv() {

        setup:
        def text = '''
            x = 1
            y = 2
            process hola {
              input:
              env x as VAR_1
              env y as 'VAR_2'

              return ''
            }
            '''

        when:
        TaskProcessor process = parse(text).run()
        EnvInParam env1 = process.taskConfig.getInputs().get(0)
        EnvInParam env2 = process.taskConfig.getInputs().get(1)

        then:
        env1.name == 'VAR_1'
        env1.inTarget == 1
        env1.inChannel.val == 1

        env2.name == 'VAR_2'
        env2.inTarget == 2
        env2.inChannel instanceof DataflowVariable
        env2.inChannel.val == 2

    }


    def tesInputEach() {

        setup:
        def text = '''

            process hola {
              input:
              each p1
              each p2 as x2
              each p3 as x3
              each p4 as x4

              return ''
            }
            '''

        def binding = [:]
        binding.p1 = 1
        binding.p2 = [1,2,3]
        binding.p3 = Nextflow.val('hola')
        binding.p4 = Nextflow.list(6,7,8)

        when:
        TaskProcessor process = parse(text, binding).run()
        EachInParam each1 = process.taskConfig.getInputs().get(0)
        EachInParam each2 = process.taskConfig.getInputs().get(1)
        EachInParam each3 = process.taskConfig.getInputs().get(2)
        EachInParam each4 = process.taskConfig.getInputs().get(3)

        then:
        each1.name == 'p1'
        each1.inTarget instanceof DataflowVariable
        each1.inChannel.val == [1]

        each2.name == 'x2'
        each2.inTarget instanceof DataflowVariable
        each2.inChannel.val == [1,2,3]

        each3.name == 'x3'
        each3.inTarget instanceof DataflowVariable
        each3.inChannel.val == ['hola']

        each4.name == 'x4'
        each4.inTarget instanceof DataflowVariable
        each4.inChannel.val == [6,7,8]
    }



    // ==============================================================
    //                  test *stdin* parameters
    // ==============================================================

    def testStdInput() {

        when:
        def text = '''
            x = 'hola'

            process hola {
              input:
              stdin x

              return ''
            }
            '''

        TaskProcessor process = parse(text).run()
        StdInParam stdin1 = process.taskConfig.getInputs().get(0)

        then:
        stdin1.name == '-'
        stdin1.inTarget == 'hola'
        stdin1.inChannel instanceof DataflowVariable
        stdin1.inChannel.val == 'hola'

    }

    def testStdInputInvalid() {

        //
        // missing 'x' variable
        //
        when:
        def text = '''

            process hola {
              input:
              stdin x

              return ''
            }
            '''
        parse(text).run()
        then:
        thrown(MissingPropertyException)

        //
        // 'as' is not supported by stdin definition
        //
        when:
        text = '''
            x = 'hola'

            process hola {
              input:
              stdin x as y
              return ''
            }
            '''
        parse(text).run()
        then:
        thrown(IllegalAccessException)

    }



    // ==============================================================
    //                  test *stdout* parameters
    // ==============================================================

    def testStdout() {

        setup:
        def text = '''
            process hola {
              output:
              stdout x
              stdout() to y

              return ''
            }
            '''

        def binding = [:]
        TaskProcessor process = parse(text, binding).run()

        when:
        StdOutParam out1 = process.taskConfig.getOutputs().get(0)
        StdOutParam out2 = process.taskConfig.getOutputs().get(1)

        // it MUST
        // - create a value out parameter named 'x'
        // - create in the script context (binding) a new variable of type DataflowQueue named 'x'
        then:
        !binding.containsKey('x')
        !binding.containsKey('y')

        out1.name == '-'
        out1.@outTarget == 'x'
        out1.outChannel instanceof DataflowQueue
        out1.outChannel == binding.x

        out2.name == '-'
        out2.@outTarget == 'y'
        out2.outChannel instanceof DataflowQueue
        out2.outChannel == binding.y

    }


    // ==============================================================
    //                  test *output* parameters
    // ==============================================================

    def testOutParam() {

        setup:
        def text = '''
            process hola {
              output:
              val x

              return ''
            }
            '''

        def binding = [:]
        TaskProcessor process = parse(text, binding).run()

        when:
        ValueOutParam out1 = process.taskConfig.getOutputs().get(0)

        // it MUST
        // - create a value out parameter named 'x'
        // - create in the script context (binding) a new variable of type DataflowQueue named 'x'
        then:
        !binding.containsKey('x')

        out1.name == 'x'
        out1.@outTarget == null
        out1.outChannel instanceof DataflowQueue
        out1.outChannel == binding.x

    }

    def testOutParam2() {

        setup:
        def text = '''
            process hola {
              output:
              val xx
              val pp to qq
              val vv to ww

              return ''
            }
            '''

        def ww = new DataflowQueue()
        def binding = [ww: ww]

        when:
        TaskProcessor process = parse(text, binding).run()
        ValueOutParam out1 = process.taskConfig.getOutputs().get(0)
        ValueOutParam out2 = process.taskConfig.getOutputs().get(1)
        ValueOutParam out3 = process.taskConfig.getOutputs().get(2)

        // it MUST

        then:
        out1.name == 'xx'
        out1.@outTarget == null
        out1.outChannel instanceof DataflowQueue
        out1.outChannel == binding.xx

        out2.name == 'pp'
        out2.@outTarget == 'qq'
        out2.outChannel instanceof DataflowQueue
        out2.outChannel == binding.qq
        !binding.containsKey('pp')

        out3.name == 'vv'
        out3.@outTarget == 'ww'
        out3.outChannel instanceof DataflowQueue
        out3.outChannel == binding.ww
        !binding.containsKey('vv')

    }


    def testOutFile1() {

        setup:
        def text = '''
            process hola {
              output:
              file blast_x

              return ''
            }
            '''

        def binding = [:]

        when:
        TaskProcessor process = parse(text, binding).run()
        FileOutParam out1 = process.taskConfig.getOutputs().get(0)

        // it MUST
        // - create a value out parameter named 'x'
        // - create in the script context (binding) a new variable of type DataflowQueue named 'x'
        then:
        !binding.containsKey('blast_x')

        out1.name == 'blast_x'
        out1.@outTarget == null
        out1.outChannel instanceof DataflowQueue
        out1.outChannel == binding.blast_x

    }


    def testOutFile2() {

        setup:
        def text = '''
            process hola {
              output:
              file x
              file y includeInputs true includeHidden true separatorChar '~' autoClose false flat true

              return ''
            }
            '''

        when:
        TaskProcessor process = parse(text).run()
        FileOutParam out1 = process.taskConfig.getOutputs().get(0)
        FileOutParam out2 = process.taskConfig.getOutputs().get(1)

        // it MUST
        // - create a value out parameter named 'x'
        // - create in the script context (binding) a new variable of type DataflowQueue named 'x'
        then:
        out1.name == 'x'
        out1.separatorChar == ':'
        out1.autoClose
        !out1.flat
        !out1.includeHidden
        !out1.includeInputs

        out2.name == 'y'
        out2.separatorChar == '~'
        out2.includeHidden
        out2.includeInputs
        !out2.autoClose
    }


    def testOutFile3() {

        setup:
        def text = '''
            process hola {
              output:
              file 'file.txt'
              file '*.fa' to p
              file '*.txt' to q

              return ''
            }
            '''

        def q = new DataflowQueue()
        def binding = [q: q]

        when:
        TaskProcessor process = parse(text, binding).run()
        FileOutParam out1 = process.taskConfig.getOutputs().get(0)
        FileOutParam out2 = process.taskConfig.getOutputs().get(1)
        FileOutParam out3 = process.taskConfig.getOutputs().get(2)

        // it MUST
        // - create a value out parameter named 'x'
        // - create in the script context (binding) a new variable of type DataflowQueue named 'file_txt'
        //   NOTE: te
        then:
        !binding.containsKey('file.txt')
        !binding.containsKey('*.txt')
        !binding.containsKey('*.fa')
        !binding.containsKey('p')
        binding.containsKey('q')

        out1.name == 'file.txt'
        out1.@outTarget == null
        out1.outChannel instanceof DataflowQueue
        out1.outChannel == binding.'file.txt'

        out2.name == '*.fa'
        out2.@outTarget == 'p'
        out2.outChannel instanceof DataflowQueue
        out2.outChannel == binding.'p'

        out3.name == '*.txt'
        out3.@outTarget == 'q'
        out3.outChannel instanceof DataflowQueue
        out3.outChannel == binding.'q'
        out3.outChannel == q

    }




    // ==============================================================
    //                  share parameters tests
    // ==============================================================

    def testSharedValueParam() {

        setup:
        def text = '''
            process hola {
              share: val x
              return ''
            }
            '''

        when:
        TaskProcessor process = parse(text).run()
        ValueSharedParam shared = process.taskConfig.getInputs().get(0)

        then:
        shared.name == 'x'
        shared.inTarget == null
        shared.@outTarget == null
        shared.inChannel instanceof DataflowVariable
        shared.inChannel.getVal() == null

    }

    def testSharedValueParam2() {

        setup:
        def text = '''
            y = 9

            process hola {
              share:
                val x
                val y

                return ''
            }
            '''

        when:
        TaskProcessor process = parse(text).run()
        ValueSharedParam shared1 = process.taskConfig.getInputs().get(0)
        ValueSharedParam shared2 = process.taskConfig.getInputs().get(1)

        then:
        process instanceof TaskProcessor
        process.taskConfig.getInputs().size() == 2

        // the 'x' shared variable binds to the default value (null)
        shared1.name == 'x'
        shared1.@inTarget == null
        shared1.@outTarget == null
        shared1.inChannel instanceof DataflowVariable
        shared1.inChannel.getVal() == null

        // the 'y' shared variable binds to the 'y' define in the script context
        shared2.name == 'y'
        shared2.@inTarget == 9
        shared2.@outTarget == null
        shared2.inChannel instanceof DataflowVariable
        shared2.inChannel.getVal() == 9

    }


    def testSharedValueParamWithUsing() {

        setup:
        def text = '''
            w = 3
            process hola {
              share:
                  val 'y' as x1
                  val { 'hola' } as x2
                  val w as x3

              return ''
            }
            '''

        when:
        TaskProcessor process = parse(text).run()
        ValueSharedParam shared1 = process.taskConfig.getInputs().get(0)
        ValueSharedParam shared2 = process.taskConfig.getInputs().get(1)
        ValueSharedParam shared3 = process.taskConfig.getInputs().get(2)

        then:
        process instanceof TaskProcessor
        process.taskConfig.getSharedDefs().size() == 3

        // 'x1' shared var binds to the value defined by 'using' keyword
        shared1.name == 'x1'
        shared1.@inTarget == 'y'
        shared1.@outTarget == null
        shared1.inChannel instanceof DataflowVariable
        shared1.inChannel.val == 'y'

        // 'x2' shared var bind to the value obtained by the closure define with the 'using' keyword
        shared2.name == 'x2'
        shared2.@inTarget instanceof Closure
        shared2.@outTarget == null
        shared2.inChannel instanceof DataflowVariable
        shared2.inChannel.val == 'hola'

        // 'x3' shared var binds 'w' variable defined in the script context
        shared3.name == 'x3'
        shared3.@inTarget == 3
        shared3.@outTarget == null
        shared3.inChannel instanceof DataflowVariable
        shared3.inChannel.val == 3

    }


    def testSharedValueParamWithInto() {

        setup:
        def text = '''
            process hola {
              share:
                  val x1 to p
                  val 22 as x2 to q
                  val 33 as x3 to w
              return ''
            }
            '''

        def w = new DataflowVariable()
        def binding = [:]
        binding.w = w

        when:
        TaskProcessor process = parse(text, binding).run()
        ValueSharedParam shared1 = process.taskConfig.getInputs().get(0)
        ValueSharedParam shared2 = process.taskConfig.getInputs().get(1)
        ValueSharedParam shared3 = process.taskConfig.getInputs().get(2)

        then:
        !binding.containsKey('p')
        !binding.containsKey('q')

        process instanceof TaskProcessor
        process.taskConfig.getSharedDefs().size() == 3

        // the 'using' declares a primitive value

        shared1.name == 'x1'
        shared1.@inTarget == null
        shared1.@outTarget == 'p'
        shared1.inChannel instanceof DataflowVariable
        shared1.outChannel instanceof DataflowVariable
        binding.containsKey('p')
        binding.get('p') == shared1.outChannel

        // the 'using' declares a closure which is resolved to the final value ('hola')
        shared2.name == 'x2'
        shared2.@inTarget == 22
        shared2.@outTarget == 'q'
        shared2.inChannel instanceof DataflowVariable

        shared2.outChannel instanceof DataflowVariable
        binding.containsKey('q')
        binding.get('q') == shared2.outChannel

        // the 'using' bind a variable declared in the script context (3)
        shared3.name == 'x3'
        shared3.@inTarget == 33
        shared3.@outTarget == 'w'
        shared3.inChannel instanceof DataflowVariable
        shared3.outChannel == w

    }


    def testSharedFile () {

        setup:
        def text = '''

            process hola {
                share:
                file 'hola'
                file x
                file y
                file z as 'file.fasta'
                file w as blast_result
                file 'str' as 'file1.txt'
                file str as 'file2.txt'

                return ''
            }
            '''

        def binding = [
                x: Paths.get('hola.txt'),
                str: 'hello world'
        ]

        when:
        TaskProcessor process = parse(text, binding).run()
        FileSharedParam shared1 = process.taskConfig.getInputs().get(0)
        FileSharedParam shared2 = process.taskConfig.getInputs().get(1)
        FileSharedParam shared3 = process.taskConfig.getInputs().get(2)
        FileSharedParam shared4 = process.taskConfig.getInputs().get(3)
        FileSharedParam shared5 = process.taskConfig.getInputs().get(4)
        FileSharedParam shared6 = process.taskConfig.getInputs().get(5)
        FileSharedParam shared7 = process.taskConfig.getInputs().get(6)

        //  define a shared file by its name
        then:
        shared1.name == null
        shared1.filePattern == null
        shared1.inChannel instanceof DataflowVariable
        shared1.inChannel.getVal() == 'hola'

        // define a shared file 'x' that binds to the same variable in the script context
        // note: the variable may hold any type not only a Path.
        shared2.name == 'x'
        shared2.filePattern == null
        shared2.inChannel.getVal() == Paths.get('hola.txt')

        // define a shared file 'y', since a variable with that name does not exist
        // binds to the default file
        shared3.name == 'y'
        shared3.filePattern == null
        shared3.inChannel.getVal() == null

        shared4.name == 'z'
        shared4.filePattern == 'file.fasta'
        shared4.inChannel.getVal() == null

        shared5.name == 'blast_result'
        shared5.filePattern == null
        shared5.inChannel.getVal() == null

        shared6.name == null
        shared6.filePattern == 'file1.txt'
        shared6.inChannel.getVal() == 'str'

        shared7.name == 'str'
        shared7.filePattern == 'file2.txt'
        shared7.inChannel.getVal() == 'hello world'

    }






}
