/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
import static test.TestParser.parseAndReturnProcess

import java.nio.file.Path

import groovyx.gpars.dataflow.DataflowQueue
import nextflow.processor.TaskContext
import nextflow.util.BlankSeparatedList
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ParamsOutTest extends Specification {

    static private createTaskContext(Map binding, Map local=[:], String name='process_name') {
        def script = new Script() {
            @Override
            Object run() { return null; }
        }

        script.setBinding(new Binding(binding))

        return new TaskContext(script, local, name)

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
              val p into q
              val 10 into t
              val 'str' into y
              val { x } into w
              val "${y}" into z

              return ''
            }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        def out0 = (ValueOutParam)process.config.getOutputs().get(0)
        def out1 = (ValueOutParam)process.config.getOutputs().get(1)
        def out2 = (ValueOutParam)process.config.getOutputs().get(2)
        def out3 = (ValueOutParam)process.config.getOutputs().get(3)
        def out4 = (ValueOutParam)process.config.getOutputs().get(4)
        def out5 = (ValueOutParam)process.config.getOutputs().get(5)

        // it MUST
        // - create a value out parameter named 'x'
        // - create in the script context (binding) a new variable of type DataflowQueue named 'x'
        then:
        process.config.getOutputs().size() == 6

        out0.name == 'x'
        out0.resolve( [x: 100] ) == 100
        out0.outChannel instanceof DataflowQueue
        binding.containsKey('x')

        out1.name == 'p'
        out1.resolve( [p: 200] ) == 200
        out1.outChannel instanceof DataflowQueue
        !binding.containsKey('p')
        binding.containsKey('q')

        out2.name == null
        out2.resolve( [:] ) == 10
        out2.outChannel instanceof DataflowQueue
        binding.containsKey('t')

        out3.name == null
        out3.resolve( [:] ) == 'str'
        out3.outChannel instanceof DataflowQueue
        binding.containsKey('y')

        out4.name == null
        out4.resolve( [x:'Hello', y:'world'] ) == 'Hello'
        out4.outChannel instanceof DataflowQueue
        binding.containsKey('w')

        out5.name == null
        out5.resolve( [x:'Hello', y:'world'] ) == 'world'
        out5.outChannel instanceof DataflowQueue
        binding.containsKey('z')

    }

    def testIntoMultipleChannels() {

        given:
        def text = '''
            process foo {
              output:
              val one into a
              val two into p, q
              file 'three' into b
              file 'four'  into x,y,z
              return ''
            }
            '''

        def binding = [:]

        when:
        def process = parseAndReturnProcess(text, binding)
        def out0 = (ValueOutParam)process.config.getOutputs().get(0)
        def out1 = (ValueOutParam)process.config.getOutputs().get(1)
        def out2 = (FileOutParam)process.config.getOutputs().get(2)
        def out3 = (FileOutParam)process.config.getOutputs().get(3)

        then:
        process.config.getOutputs().size() == 4

        out0.name == 'one'
        out0.getOutChannels().size()==1
        out0.getOutChannels().get(0) instanceof DataflowQueue

        out1.name == 'two'
        out1.getOutChannels().size()==2
        out1.getOutChannels().get(0) instanceof DataflowQueue
        out1.getOutChannels().get(1) instanceof DataflowQueue

        out2.name == null
        out2.getOutChannels().size()==1
        out2.getOutChannels().get(0) instanceof DataflowQueue

        out3.name == null
        out3.getOutChannels().size()==3
        out3.getOutChannels().get(0) instanceof DataflowQueue
        out3.getOutChannels().get(1) instanceof DataflowQueue
        out3.getOutChannels().get(2) instanceof DataflowQueue
    }

    def testFileOutParams() {

        setup:
        def text = '''
            process hola {
              output:
              file x
              file 'y' mode flatten
              file p into q mode standard

              return ''
            }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        def out1 = process.config.getOutputs().get(0)
        def out2 = process.config.getOutputs().get(1)
        def out3 = process.config.getOutputs().get(2)

        // it MUST
        // - create a value out parameter named 'x'
        // - create in the script context (binding) a new variable of type DataflowQueue named 'x'
        then:
        process.config.getOutputs().size() == 3

        out1.class == FileOutParam
        out1.name == 'x'
        out1.outChannel instanceof DataflowQueue
        out1.outChannel == binding.'x'
        out1.mode == BasicMode.standard

        out2.class == FileOutParam
        out2.name == null
        out2.outChannel == null
        out2.mode == BasicMode.flatten

        out3.class == FileOutParam
        out3.name == 'p'
        out3.outChannel instanceof DataflowQueue
        out3.outChannel == binding.q
        out3.mode == BasicMode.standard
        !binding.containsKey('p')
    }

    def testFileOutParamsWithVariables() {

        setup:
        def text = '''

            process hola {
              output:
              file "${x}_name" into channel1
              file "${x}_${y}.fa" into channel2
              file "simple.txt" into channel3
              file "${z}.txt:sub/dir/${x}.fa" into channel4
              set "${z}.txt:${x}.fa" into channel5
              set file("${z}.txt:${x}.fa") into channel6

              return ''
            }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)
        def ctx = [x: 'hola', y:99, z:'script_file']

        when:
        FileOutParam out0 = process.config.getOutputs().get(0)
        FileOutParam out1 = process.config.getOutputs().get(1)
        FileOutParam out2 = process.config.getOutputs().get(2)
        FileOutParam out3 = process.config.getOutputs().get(3)
        SetOutParam out4 = process.config.getOutputs().get(4)
        SetOutParam out5 = process.config.getOutputs().get(5)

        then:
        process.config.getOutputs().size() == 6

        out0.name == null
        out0.getFilePatterns(ctx,null) == ['hola_name']
        out0.outChannel instanceof DataflowQueue
        out0.outChannel == binding.channel1
        out0.isDynamic()

        out1.name == null
        out1.getFilePatterns(ctx,null) == ['hola_99.fa']
        out1.outChannel instanceof DataflowQueue
        out1.outChannel == binding.channel2
        out1.isDynamic()

        out2.name == null
        out2.getFilePatterns(ctx,null) == ['simple.txt']
        out2.outChannel instanceof DataflowQueue
        out2.outChannel == binding.channel3
        !out2.isDynamic()

        out3.name == null
        out3.getFilePatterns(ctx,null) == ['script_file.txt','sub/dir/hola.fa']
        out3.outChannel instanceof DataflowQueue
        out3.outChannel == binding.channel4
        out3.isDynamic()

        out4.name == 'setoutparam<4>'
        out4.outChannel instanceof DataflowQueue
        out4.outChannel == binding.channel5
        out4.inner[0] instanceof FileOutParam
        (out4.inner[0] as FileOutParam) .getFilePatterns(ctx,null) == ['script_file.txt','hola.fa']
        (out4.inner[0] as FileOutParam) .isDynamic()

        out5.name == 'setoutparam<5>'
        out5.outChannel instanceof DataflowQueue
        out5.outChannel == binding.channel6
        out5.inner[0] instanceof FileOutParam
        (out5.inner[0] as FileOutParam) .getFilePatterns(ctx,null) == ['script_file.txt','hola.fa']
        (out5.inner[0] as FileOutParam) .isDynamic()

    }

    def testFileOutFileCollection () {

        given:
        def text = '''

            process hola {
              output:
              file(x) into channel1
              set file(y) into channel2

              return ''
            }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        def list_x = new BlankSeparatedList(['a.txt' as Path, 'b.txt' as Path, 'c.txt' as Path])
        def list_y = ['one.txt' as Path, 'two.txt' as Path, 'three.txt' as Path]
        def ctx = [x: list_x, y: list_y]

        when:
        FileOutParam out0 = process.config.getOutputs().get(0)
        SetOutParam out1 = process.config.getOutputs().get(1)

        then:
        out0.getFilePatterns(ctx,null) == ['a.txt', 'b.txt', 'c.txt' ]
        out1.inner[0].getFilePatterns(ctx,null) == ['one.txt', 'two.txt', 'three.txt']

    }

    def testFileOutWithGString2 () {

        setup:
        def text = '''

            process hola {
              output:
              file x
              file "$y" into q
              set file(z) into p
              file u
              file "$v"
              file 'w'

              return ''
            }
            '''

        def binding = [ x: 'hola', y:'hola_2', z: 'hola_z', v:'file_v' ]
        def process = parseAndReturnProcess(text, binding)

        when:
        FileOutParam out0 = process.config.getOutputs().get(0)
        FileOutParam out1 = process.config.getOutputs().get(1)
        SetOutParam out2 = process.config.getOutputs().get(2)
        FileOutParam out3 = process.config.getOutputs().get(3)
        FileOutParam out4 = process.config.getOutputs().get(4)
        FileOutParam out5 = process.config.getOutputs().get(5)

        then:
        out0.name == 'x'
        out0.getFilePatterns(binding,null) == ['hola']
        out0.getOutChannels().get(0) instanceof DataflowQueue
        out0.getOutChannels().get(0) == binding.x

        out1.name == null
        out1.getFilePatterns(binding,null) == ['hola_2']
        out1.getOutChannels().get(0) instanceof DataflowQueue
        out1.getOutChannels().get(0) == binding.q

        out2.inner[0] instanceof FileOutParam
        (out2.inner[0] as FileOutParam).name == 'z'
        (out2.inner[0] as FileOutParam).getFilePatterns(binding,null) == ['hola_z']

        out3.name == 'u'
        out3.getFilePatterns(binding,null) == ['u']
        out3.getOutChannels().get(0) instanceof DataflowQueue
        out3.getOutChannels().get(0) == binding.u

        out4.name == null
        out4.getFilePatterns(binding,null) == ['file_v']
        out4.getOutChannels().size()==0

        out5.name == null
        out5.getFilePatterns(binding,null) == ['w']
        out4.getOutChannels().size()==0
    }


    def testFileOutWithClosure() {

        setup:
        def text = '''

            process hola {
              output:
              file { "${x}_name" } into channel1
              file { "${params.fileName}_${y}.fa" } into channel2
              set file({ "${z}.txt" }) into channel3

              return ''
            }
            '''

        def binding = [x: 'hola', y:99, z:'script_file', params: [fileName: 'hello']]
        def process = parseAndReturnProcess(text, binding)

        when:
        FileOutParam out1 = process.config.getOutputs().get(0)
        FileOutParam out2 = process.config.getOutputs().get(1)
        SetOutParam out3 = process.config.getOutputs().get(2)

        then:
        out1.name == null
        out1.getFilePatterns(binding,null) == ['hola_name']
        out1.outChannel instanceof DataflowQueue
        out1.outChannel == binding.channel1
        out1.isDynamic()

        out2.name == null
        out2.getFilePatterns(binding,null) == ['hello_99.fa']
        out2.outChannel instanceof DataflowQueue
        out2.outChannel == binding.channel2
        out2.isDynamic()

        out3.outChannel instanceof DataflowQueue
        out3.outChannel == binding.channel3
        out3.inner[0] instanceof FileOutParam
        (out3.inner[0] as FileOutParam) .getFilePatterns(binding,null) == ['script_file.txt']
        (out3.inner[0] as FileOutParam) .isDynamic()

    }

    def testFileOutWithParams() {

        setup:
        def text = '''

            process hola {
              output:
              file x into ch

              file x maxDepth 5 into ch
              file x hidden true into ch
              file x followLinks false into ch
              file x type 'file' into ch
              file x separatorChar '#' into ch

              file x hidden false into ch
              file x followLinks true into ch
              file x type 'dir' into ch
              file x glob false into ch
              file x optional true into ch

              return ''
            }
            '''

        def process = parseAndReturnProcess(text, [:])

        when:
        FileOutParam out0 = process.config.getOutputs().get(0)
        FileOutParam out1 = process.config.getOutputs().get(1)
        FileOutParam out2 = process.config.getOutputs().get(2)
        FileOutParam out3 = process.config.getOutputs().get(3)
        FileOutParam out4 = process.config.getOutputs().get(4)
        FileOutParam out5 = process.config.getOutputs().get(5)
        FileOutParam out6 = process.config.getOutputs().get(6)
        FileOutParam out7 = process.config.getOutputs().get(7)
        FileOutParam out8 = process.config.getOutputs().get(8)
        FileOutParam out9 = process.config.getOutputs().get(9)
        FileOutParam out10 = process.config.getOutputs().get(10)

        then:
        out0.maxDepth == null
        !out0.hidden
        out0.followLinks
        out0.type == null
        out0.separatorChar == ':'
        out0.glob
        !out0.optional

        out1.maxDepth == 5
        out2.hidden
        !out3.followLinks
        out4.type == 'file'
        out5.separatorChar == '#'
        !out6.hidden
        out7.followLinks
        out8.type == 'dir'
        out9.glob == false
        out10.optional

    }

    def testFileWithoutInto() {
        setup:
        def text = '''
            process hola {
              output:
                file 'x'
                file y
                file 'z' into channel_z

              return ''
            }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        FileOutParam out0 = process.config.getOutputs().get(0)
        FileOutParam out1 = process.config.getOutputs().get(1)
        FileOutParam out2 = process.config.getOutputs().get(2)

        then:
        out0.name == null
        out0.getFilePatterns(binding,null) == ['x']
        out0.getOutChannels().size()==0

        out1.name == 'y'
        out1.getFilePatterns(binding,null) == ['y']
        out1.getOutChannels().get(0) instanceof DataflowQueue
        out1.getOutChannels().get(0) == binding.y

        out2.name == null
        out2.getFilePatterns(binding,null) == ['z']
        out2.getOutChannels().get(0) instanceof DataflowQueue
        out2.getOutChannels().get(0) == binding.channel_z
    }


    def testSetOutParams() {

        setup:
        def text = '''
            process hola {
              output:
                set(x) into p
                set(y, '-', '*.fa') into q mode flatten
                set(stdout, z) into t mode combine

              return ''
            }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        SetOutParam out1 = process.config.getOutputs().get(0)
        SetOutParam out2 = process.config.getOutputs().get(1)
        SetOutParam out3 = process.config.getOutputs().get(2)

        then:
        process.config.getOutputs().size() == 3

        out1.outChannel instanceof DataflowQueue
        out1.outChannel == binding.p
        out1.inner.size() == 1
        out1.inner[0] instanceof ValueOutParam
        out1.inner[0].name == 'x'
        out1.inner[0].index == 0
        out1.mode == BasicMode.standard

        out2.outChannel instanceof DataflowQueue
        out2.outChannel == binding.q
        out2.inner[0] instanceof ValueOutParam
        out2.inner[0].name == 'y'
        out2.inner[0].index == 1
        out2.inner[1] instanceof StdOutParam
        out2.inner[1].name == '-'
        out2.inner[1].index == 1
        out2.inner[2] instanceof FileOutParam
        out2.inner[2].name == null
        out2.inner[2].filePattern == '*.fa'
        out2.inner[2].index == 1
        out2.inner.size() ==3
        out2.mode == BasicMode.flatten

        out3.outChannel instanceof DataflowQueue
        out3.outChannel == binding.t
        out3.inner.size() == 2
        out3.inner[0] instanceof StdOutParam
        out3.inner[0].name == '-'
        out3.inner[0].index == 2
        out3.inner[1] instanceof ValueOutParam
        out3.inner[1].name == 'z'
        out3.inner[1].index == 2
        out3.mode == SetOutParam.CombineMode.combine

    }

    def testSetOutParams2() {

        setup:
        def text = '''
            process hola {
              output:
                set val(x) into p
                set val(y), stdout, file('*.fa') into q mode flatten
                set stdout, val(z) into t mode combine

              return ''
            }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        SetOutParam out0 = process.config.getOutputs().get(0)
        SetOutParam out1 = process.config.getOutputs().get(1)
        SetOutParam out2 = process.config.getOutputs().get(2)

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
        out2.mode == SetOutParam.CombineMode.combine

    }

    def testSetOutValues() {

        setup:
        def text = '''
            process hola {
              output:
                set val(x), val('x') into p
                set val(1), val('2') into q
                set val("$foo"), val { bar } into t

              return ''
            }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        SetOutParam out0 = process.config.getOutputs().get(0)
        SetOutParam out1 = process.config.getOutputs().get(1)
        SetOutParam out2 = process.config.getOutputs().get(2)

        then:
        process.config.getOutputs().size() == 3

        // first set
        out0.outChannel instanceof DataflowQueue
        out0.outChannel == binding.p
        out0.inner.size() == 2
        // val(x)
        out0.inner[0] instanceof ValueOutParam
        out0.inner[0].name == 'x'
        out0.inner[0].index == 0
        out0.inner[0].resolve([x: 'hello']) == 'hello'
        // val('x')
        out0.inner[1] instanceof ValueOutParam
        out0.inner[1].name == null
        out0.inner[1].index == 0
        out0.inner[1].resolve([:]) == 'x'

        // -- second set
        out1.outChannel instanceof DataflowQueue
        out1.outChannel == binding.q
        out1.inner.size() == 2
        // val(1)
        out1.inner[0] instanceof ValueOutParam
        out1.inner[0].name == null
        out1.inner[0].index == 1
        out1.inner[0].resolve([:]) == 1
        // val('2')
        out1.inner[1] instanceof ValueOutParam
        out1.inner[1].name == null
        out1.inner[1].index == 1
        out1.inner[1].resolve([:]) == '2'

        // -- third set
        out2.outChannel instanceof DataflowQueue
        out2.outChannel == binding.t
        out2.inner.size() == 2
        // val("$foo")
        out2.inner[0] instanceof ValueOutParam
        out2.inner[0].name == null
        out2.inner[0].index == 2
        out2.inner[0].resolve([foo:'hello', bar:'world']) == 'hello'
        // val { bar }
        out2.inner[1] instanceof ValueOutParam
        out2.inner[1].name == null
        out2.inner[1].index == 2
        out2.inner[1].resolve([foo:'hello', bar:'world']) == 'world'

    }

    def testSetOutWithoutInto() {

        setup:
        def text = '''
            process hola {
              output:
              set val(X), file('Y')

              return ''
            }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        SetOutParam out0 = process.config.getOutputs().get(0)

        then:
        process.config.getOutputs().size() == 1

        // first set
        out0.getOutChannels().size()==0

        out0.inner[0] instanceof ValueOutParam
        out0.inner[0].name == 'X'
        out0.inner[0].index == 0
        out0.inner[0].mapIndex == 0

        out0.inner[1] instanceof FileOutParam
        out0.inner[1].name == null
        out0.inner[1].filePattern == 'Y'
        out0.inner[1].index == 0
        out0.inner[1].mapIndex == 1

    }

    def testStdOut() {

        setup:
        def text = '''
            process hola {
              output:
              stdout into p
              stdout into (q)
              stdout into (x,y,z)

              return ''
            }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        def out0 = (StdOutParam)process.config.getOutputs().get(0)
        def out1 = (StdOutParam)process.config.getOutputs().get(1)
        def out2 = (StdOutParam)process.config.getOutputs().get(2)

        then:
        process.config.getOutputs().size() == 3

        out0.getOutChannels()[0].is binding.p
        out1.getOutChannels()[0].is binding.q
        out2.getOutChannels()[0].is binding.x
        out2.getOutChannels()[1].is binding.y
        out2.getOutChannels()[2].is binding.z

    }


    def testOutList () {

        setup:
        def bind = new Binding()
        def outs = new OutputsList()

        when:
        def v1 = new ValueOutParam(bind,outs)
        def o1 = new StdOutParam(bind,outs)
        def v2 = new ValueOutParam(bind,outs)
        def s1 = new SetOutParam(bind,outs)
        def s2 = new SetOutParam(bind,outs)

        then:
        outs.size() == 5
        outs.ofType(StdOutParam) == [o1]
        outs.ofType(ValueOutParam) == [v1,v2]
        outs.ofType(SetOutParam,StdOutParam) == [o1,s1,s2]

    }


    def testModeParam() {

        setup:
        def p = new SetOutParam(new Binding(), [])
        when:
        p.mode(value)
        then:
        p.getMode() == expected

        where:
        value                       | expected
        'combine'                   | SetOutParam.CombineMode.combine
        new TokenVar('combine')     | SetOutParam.CombineMode.combine
        'flatten'                   | BasicMode.flatten
        new TokenVar('flatten')     | BasicMode.flatten

    }

    def testWrongMode() {

        when:
        def p = new SetOutParam(new Binding(), [])
        p.mode('unknown')
        then:
        thrown(IllegalArgumentException)

    }

    def testDefaultMode() {

        setup:
        def bind = new Binding()
        def list = []

        expect:
        new StdOutParam(bind, list).mode == BasicMode.standard
        new ValueOutParam(bind, list).mode == BasicMode.standard
        new FileOutParam(bind, list).mode == BasicMode.standard
        new SetOutParam(bind, list).mode == BasicMode.standard

    }

    def 'should fetch output value' () {

        given:
        def param = new ValueOutParam(new Binding(), [])

        /*
         * test input:
         * val x
         */
        when:
        param.target = new TokenVar('x')
        then:
        param.resolve(createTaskContext([x:'foo'])) == 'foo'

        when:
        param.resolve(createTaskContext([z:'foo']))
        then:
        def error = thrown(MissingPropertyException)
        error.property == 'x'

        /*
         * test input:
         * val { str }
         */
        when:
        param.target = { str }
        then:
        param.resolve( createTaskContext([str:'foo']) ) == 'foo'

        when:
        param.resolve(createTaskContext([:]))
        then:
        error = thrown(MissingPropertyException)
        error.property == 'str'

        /*
         * test input:
         * val "${params.data}"
         */
        when:
        param.target = "${->params.data}"
        then:
        param.resolve( createTaskContext([params:[data:99]]) ) == '99'

        when:
        param.resolve(createTaskContext([:]))
        then:
        error = thrown(MissingPropertyException)
        error.property == 'params'
    }

}
