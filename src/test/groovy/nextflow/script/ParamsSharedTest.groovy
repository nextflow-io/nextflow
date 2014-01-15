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

package nextflow.script
import static test.TestParser.parse

import groovyx.gpars.dataflow.DataflowVariable
import nextflow.processor.TaskProcessor
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ParamsSharedTest extends Specification {


    def testSharedValue() {


        setup:
        def text = '''

            process hola {
              share:
              val x
              val y into z
              val q into q

              return ''
            }
            '''

        def binding = [:]
        when:
        TaskProcessor process = parse(text, binding).run()
        ValueSharedParam share1 = process.taskConfig.getInputs().get(0)
        ValueSharedParam share2 = process.taskConfig.getInputs().get(1)
        ValueSharedParam share3 = process.taskConfig.getInputs().get(2)

        then:
        process.taskConfig.getInputs().size() == 3

        share1.name == 'x'
        share1.inChannel instanceof DataflowVariable
        share1.inChannel.val == null
        share1.outChannel == null
        !binding.containsKey('x')

        share2.name == 'y'
        share2.inChannel instanceof DataflowVariable
        share2.inChannel.val == null
        share2.outChannel instanceof DataflowVariable
        binding.z == share2.outChannel
        !binding.containsKey('y')

        share3.name == 'q'
        share3.inChannel instanceof DataflowVariable
        share3.inChannel.val == null
        share3.outChannel instanceof DataflowVariable
        binding.containsKey('q')
        binding.q == share3.outChannel

    }

    def testSharedTargetChannelNameCollision() {

        setup:
        def text = '''
            x = 1

            process hola {
              share:
              val x into x

              return ''
            }
            '''

        def binding = [:]
        when:
        TaskProcessor process = parse(text, binding).run()
        process.taskConfig.getInputs().get(0).inChannel

        then:
        thrown(IllegalArgumentException)

    }

    def testSharedFromValue() {


        setup:
        def text = '''
            x = 1
            y = 2
            v = 3

            process hola {
              share:
              val x
              val y into z
              val w from v into q

              return ''
            }
            '''

        def binding = [:]
        when:
        TaskProcessor process = parse(text, binding).run()
        ValueSharedParam share1 = process.taskConfig.getInputs().get(0)
        ValueSharedParam share2 = process.taskConfig.getInputs().get(1)
        ValueSharedParam share3 = process.taskConfig.getInputs().get(2)

        then:
        process.taskConfig.getInputs().size() == 3

        share1.name == 'x'
        share1.inChannel instanceof DataflowVariable
        share1.inChannel.val == 1
        share1.outChannel == null

        share2.name == 'y'
        share2.inChannel instanceof DataflowVariable
        share2.inChannel.val == 2
        share2.outChannel instanceof DataflowVariable
        binding.get('z') == share2.outChannel

        share3.name == 'w'
        share3.inChannel instanceof DataflowVariable
        share3.inChannel.val == 3
        share3.outChannel instanceof DataflowVariable
        binding.get('q') == share3.outChannel
        !binding.containsKey('w')

    }



    def testSharedFromValueLiterals() {

        setup:
        def text = '''
            process hola {

              share:
              val x from 1
              val y from { 'str' }

              return ''
            }
            '''

        when:
        TaskProcessor process = parse(text).run()
        ValueSharedParam share1 = process.taskConfig.getInputs().get(0)
        ValueSharedParam share2 = process.taskConfig.getInputs().get(1)

        then:
        process.taskConfig.getInputs().size() ==2

        share1.name == 'x'
        share1.inChannel.val == 1
        share1.outChannel == null

        share2.name == 'y'
        share2.inChannel.val == 'str'
        share2.outChannel == null

    }

    def testSharedFromFile() {


        setup:
        def text = '''
            x = 1
            p = 2

            process hola {
              share:
              file alpha
              file x
              file q from p
              file 'alpha.txt' from p
              file q:'beta.txt' from p
              file 'hello.txt' into z from p

              return ''
            }
            '''

        def binding = [:]
        when:
        TaskProcessor process = parse(text, binding).run()
        FileSharedParam share0 = process.taskConfig.getInputs().get(0)
        FileSharedParam share1 = process.taskConfig.getInputs().get(1)
        FileSharedParam share2 = process.taskConfig.getInputs().get(2)
        FileSharedParam share3 = process.taskConfig.getInputs().get(3)
        FileSharedParam share4 = process.taskConfig.getInputs().get(4)
        FileSharedParam share5 = process.taskConfig.getInputs().get(5)

        then:
        process.taskConfig.getInputs().size() == 6

        share0.name == 'alpha'
        share0.inChannel instanceof DataflowVariable
        share0.inChannel.val == null
        share0.outChannel == null

        share1.name == 'x'
        share1.inChannel instanceof DataflowVariable
        share1.inChannel.val == 1    // bind automatically to 'x' even tough 'from' is omitter
        share1.outChannel == null

        share2.name == 'q'
        share2.filePattern == '*'
        share2.inChannel instanceof DataflowVariable
        share2.inChannel.val == 2
        share2.outChannel == null

        share3.name == 'alpha.txt'
        share3.filePattern == 'alpha.txt'
        share3.inChannel instanceof DataflowVariable
        share3.inChannel.val == 2
        share2.outChannel == null

        share4.name == 'q'
        share4.filePattern == 'beta.txt'
        share4.inChannel instanceof DataflowVariable
        share4.inChannel.val == 2
        share4.outChannel == null

        share5.name == 'hello.txt'
        share5.filePattern == 'hello.txt'
        share5.inChannel instanceof DataflowVariable
        share5.inChannel.val == 2
        share5.outChannel instanceof DataflowVariable
        binding.get('z') == share5.outChannel


    }


    def testDefaultMode() {

        setup:
        def bind = new Binding()
        def list = []

        expect:
        new FileSharedParam(bind, list).mode == BasicMode.standard
        new ValueSharedParam(bind, list).mode == BasicMode.standard
    }

    def testModeParam() {

        setup:
        def p = new ValueSharedParam(new Binding(), [])
        when:
        p.mode(value)
        then:
        p.getMode() == expected

        where:
        value                       | expected
        'flatten'                   | BasicMode.flatten
        new ScriptVar('flatten')    | BasicMode.flatten
        'standard'                  | BasicMode.standard
        new ScriptVar('standard')   | BasicMode.standard

    }

    def testWrongMode() {

        when:
        def p = new ValueSharedParam(new Binding(), [])
        p.mode('unknown')
        then:
        thrown(IllegalArgumentException)

    }



}
