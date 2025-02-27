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

import static test.TestParser.*

import java.nio.file.Path

import groovyx.gpars.dataflow.DataflowVariable
import nextflow.processor.TaskContext
import nextflow.script.TokenVar
import nextflow.util.BlankSeparatedList
import test.Dsl2Spec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ParamsOutTest extends Dsl2Spec {

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
              val p
              val 10
              val 'str' 
              val { x }
              val "${y}"
              val x.y 

              return ''
            }
            
            workflow {
              hola()
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
        def out6 = (ValueOutParam)process.config.getOutputs().get(6)

        // it MUST
        // - create a value out parameter named 'x'
        // - create in the script context (binding) a new variable of type DataflowQueue named 'x'
        then:
        process.config.getOutputs().size() == 7

        out0.name == 'x'
        out0.resolve( [x: 100] ) == 100
        out0.outChannel instanceof DataflowVariable

        out1.name == 'p'
        out1.resolve( [p: 200] ) == 200
        out1.outChannel instanceof DataflowVariable

        out2.name == null
        out2.resolve( [:] ) == 10
        out2.outChannel instanceof DataflowVariable

        out3.name == null
        out3.resolve( [:] ) == 'str'
        out3.outChannel instanceof DataflowVariable

        out4.name == null
        out4.resolve( [x:'Hello', y:'world'] ) == 'Hello'
        out4.outChannel instanceof DataflowVariable

        out5.name == null
        out5.resolve( [x:'Hello', y:'world'] ) == 'world'
        out5.outChannel instanceof DataflowVariable

        out6.name == null
        out6.resolve( [x:[y:'Ciao']] ) == 'Ciao'
        out6.outChannel instanceof DataflowVariable
    }

    def testIntoMultipleChannels() {

        given:
        def text = '''
            process foo {
              output:
              val one 
              file 'two'
              return ''
            }
            
            workflow { foo() }
            '''

        def binding = [:]

        when:
        def process = parseAndReturnProcess(text, binding)
        def out0 = (ValueOutParam)process.config.getOutputs().get(0)
        def out1 = (FileOutParam)process.config.getOutputs().get(1)

        then:
        process.config.getOutputs().size() == 2

        out0.name == 'one'
        out0.getOutChannel() instanceof DataflowVariable

        out1.name == null
        out1.getOutChannel() instanceof DataflowVariable

    }

    def testFileOutParams() {

        setup:
        def text = '''
            process hola {
              output:
              file x
              file 'y'
              file p

              return ''
            }
            
            workflow { hola() }
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
        out1.outChannel instanceof DataflowVariable

        out2.class == FileOutParam
        out2.name == null
        out2.outChannel instanceof DataflowVariable

        out3.class == FileOutParam
        out3.name == 'p'
        out3.outChannel instanceof DataflowVariable
    }

    def testFileOutParamsWithVariables() {

        setup:
        def text = '''

            process hola {
              output:
              file "${x}_name" 
              file "${x}_${y}.fa" 
              file "simple.txt" 
              file "${z}.txt:sub/dir/${x}.fa" 
              tuple file("${z}.txt:${x}.fa") 
              tuple path("${z}.txt:${x}.fa") 
              file meta.id  
              file "$meta.id" 
              return ''
            }
            
            workflow { hola() }
            '''

        def binding = [x: 'hola', y:99, z:'script_file', meta: [id:'hello.txt']]
        def process = parseAndReturnProcess(text, binding)
        def ctx = binding
        
        when:
        FileOutParam out0 = process.config.getOutputs().get(0)
        FileOutParam out1 = process.config.getOutputs().get(1)
        FileOutParam out2 = process.config.getOutputs().get(2)
        FileOutParam out3 = process.config.getOutputs().get(3)
        TupleOutParam out4 = process.config.getOutputs().get(4)
        TupleOutParam out5 = process.config.getOutputs().get(5)
        FileOutParam out6 = process.config.getOutputs().get(6)
        FileOutParam out7 = process.config.getOutputs().get(7)

        then:
        process.config.getOutputs().size() == 8

        out0.name == null
        out0.getFilePatterns(ctx,null) == ['hola_name']
        out0.outChannel instanceof DataflowVariable
        out0.isDynamic()

        out1.name == null
        out1.getFilePatterns(ctx,null) == ['hola_99.fa']
        out1.outChannel instanceof DataflowVariable
        out1.isDynamic()

        out2.name == null
        out2.getFilePatterns(ctx,null) == ['simple.txt']
        out2.outChannel instanceof DataflowVariable
        !out2.isDynamic()

        out3.name == null
        out3.getFilePatterns(ctx,null) == ['script_file.txt','sub/dir/hola.fa']
        out3.outChannel instanceof DataflowVariable
        out3.isDynamic()

        out4.name == 'tupleoutparam<4>'
        out4.outChannel instanceof DataflowVariable
        out4.inner[0] instanceof FileOutParam
        (out4.inner[0] as FileOutParam) .getFilePatterns(ctx,null) == ['script_file.txt','hola.fa']
        (out4.inner[0] as FileOutParam) .isDynamic()

        out5.name == 'tupleoutparam<5>'
        out5.outChannel instanceof DataflowVariable
        out5.inner[0] instanceof FileOutParam
        (out5.inner[0] as FileOutParam) .getFilePatterns(ctx,null) == ['script_file.txt:hola.fa']
        (out5.inner[0] as FileOutParam) .isDynamic()

        out6.name == null
        out6.getFilePatterns(ctx,null) == ['hello.txt']
        out6.outChannel instanceof DataflowVariable
        out6.isDynamic()

        out7.name == null
        out7.getFilePatterns(ctx,null) == ['hello.txt']
        out7.outChannel instanceof DataflowVariable
        out7.isDynamic()
    }

    def testFileOutFileCollection () {

        given:
        def text = '''

            process hola {
              output:
              file(x)
              tuple file(y)

              return ''
            }
            
            workflow { hola() }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        def list_x = new BlankSeparatedList(['a.txt' as Path, 'b.txt' as Path, 'c.txt' as Path])
        def list_y = ['one.txt' as Path, 'two.txt' as Path, 'three.txt' as Path]
        def ctx = [x: list_x, y: list_y]

        when:
        FileOutParam out0 = process.config.getOutputs().get(0)
        TupleOutParam out1 = process.config.getOutputs().get(1)

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
              file "$y"
              tuple file(z)
              file u
              file "$v"
              file 'w'

              return ''
            }
            
            workflow { hola() }
            '''

        def binding = [ x: 'hola', y:'hola_2', z: 'hola_z', v:'file_v' ]
        def process = parseAndReturnProcess(text, binding)

        when:
        FileOutParam out0 = process.config.getOutputs().get(0)
        FileOutParam out1 = process.config.getOutputs().get(1)
        TupleOutParam out2 = process.config.getOutputs().get(2)
        FileOutParam out3 = process.config.getOutputs().get(3)
        FileOutParam out4 = process.config.getOutputs().get(4)
        FileOutParam out5 = process.config.getOutputs().get(5)

        then:
        out0.name == 'x'
        out0.getFilePatterns(binding,null) == ['hola']
        out0.getOutChannel() instanceof DataflowVariable

        out1.name == null
        out1.getFilePatterns(binding,null) == ['hola_2']
        out1.getOutChannel() instanceof DataflowVariable

        out2.inner[0] instanceof FileOutParam
        (out2.inner[0] as FileOutParam).name == 'z'
        (out2.inner[0] as FileOutParam).getFilePatterns(binding,null) == ['hola_z']

        out3.name == 'u'
        out3.getFilePatterns(binding,null) == ['u']
        out3.getOutChannel() instanceof DataflowVariable

        out4.name == null
        out4.getFilePatterns(binding,null) == ['file_v']
        out4.getOutChannel() instanceof DataflowVariable

        out5.name == null
        out5.getFilePatterns(binding,null) == ['w']
        out4.getOutChannel() instanceof DataflowVariable
    }


    def testFileOutWithClosure() {

        setup:
        def text = '''

            process hola {
              output:
              file { "${x}_name" } 
              file { "${params.fileName}_${y}.fa" } 
              tuple file({ "${z}.txt" }) 

              return ''
            }
            
            workflow { hola() }
            '''

        def binding = [x: 'hola', y:99, z:'script_file', params: [fileName: 'hello']]
        def process = parseAndReturnProcess(text, binding)

        when:
        FileOutParam out1 = process.config.getOutputs().get(0)
        FileOutParam out2 = process.config.getOutputs().get(1)
        TupleOutParam out3 = process.config.getOutputs().get(2)

        then:
        out1.name == null
        out1.getFilePatterns(binding,null) == ['hola_name']
        out1.outChannel instanceof DataflowVariable
        out1.isDynamic()

        out2.name == null
        out2.getFilePatterns(binding,null) == ['hello_99.fa']
        out2.outChannel instanceof DataflowVariable
        out2.isDynamic()

        out3.outChannel instanceof DataflowVariable
        out3.inner[0] instanceof FileOutParam
        (out3.inner[0] as FileOutParam) .getFilePatterns(binding,null) == ['script_file.txt']
        (out3.inner[0] as FileOutParam) .isDynamic()

    }

    def testFileOutWithParams() {

        setup:
        def text = '''

            process hola {
              output:
              file x 

              path x, maxDepth: 5
              path x, hidden: true 
              path x, followLinks: false 
              path x, type: 'file' 
              path x, separatorChar: '#' 

              path x, hidden: false 
              path x, followLinks: true 
              path x, type: 'dir' 
              path x, glob: false 
              path x, optional: true 

              return ''
            }
            
            workflow { hola() }
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

    def testTupleOutParams() {

        setup:
        def text = '''
            process hola {
              output:
                tuple val(x)
                tuple val(y), stdout, file('*.fa') 
                tuple stdout, val(z)

              return ''
            }
            
            workflow { hola() }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        TupleOutParam out1 = process.config.getOutputs().get(0)
        TupleOutParam out2 = process.config.getOutputs().get(1)
        TupleOutParam out3 = process.config.getOutputs().get(2)

        then:
        process.config.getOutputs().size() == 3

        out1.outChannel instanceof DataflowVariable
        out1.inner.size() == 1
        out1.inner[0] instanceof ValueOutParam
        out1.inner[0].name == 'x'
        out1.inner[0].index == 0

        out2.outChannel instanceof DataflowVariable
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

        out3.outChannel instanceof DataflowVariable
        out3.inner.size() == 2
        out3.inner[0] instanceof StdOutParam
        out3.inner[0].name == '-'
        out3.inner[0].index == 2
        out3.inner[1] instanceof ValueOutParam
        out3.inner[1].name == 'z'
        out3.inner[1].index == 2

    }

    def testTupleOutValues() {

        setup:
        def text = '''
            process hola {
              output:
                tuple val(x), val('x')
                tuple val(1), val('2')
                tuple val("$foo"), val { bar }

              return ''
            }
            
            workflow { hola() }  
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        TupleOutParam out0 = process.config.getOutputs().get(0)
        TupleOutParam out1 = process.config.getOutputs().get(1)
        TupleOutParam out2 = process.config.getOutputs().get(2)

        then:
        process.config.getOutputs().size() == 3

        // first tuple
        out0.outChannel instanceof DataflowVariable
        out0.inner.size() == 2
        // val(x)
        out0.inner[0] instanceof ValueOutParam
        out0.inner[0].index == 0
        out0.inner[0].resolve([x: 'hello']) == 'hello'
        // val('x')
        out0.inner[1] instanceof ValueOutParam
        out0.inner[1].name == null
        out0.inner[1].index == 0
        out0.inner[1].resolve([:]) == 'x'

        // -- second tuple
        out1.outChannel instanceof DataflowVariable
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

        // -- third tuple
        out2.outChannel instanceof DataflowVariable
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

    def testStdOut() {

        setup:
        def text = '''
            process hola {
              output:
              stdout

              return ''
            }
            
            workflow { hola() }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        def out0 = (StdOutParam)process.config.getOutputs().get(0)

        then:
        process.config.getOutputs().size() == 1

        out0.getOutChannel() instanceof DataflowVariable

    }


    def testOutList () {

        setup:
        def bind = new Binding()
        def outs = new OutputsList()

        when:
        def v1 = new ValueOutParam(bind,outs)
        def o1 = new StdOutParam(bind,outs)
        def v2 = new ValueOutParam(bind,outs)
        def s1 = new TupleOutParam(bind,outs)
        def s2 = new TupleOutParam(bind,outs)

        then:
        outs.size() == 5
        outs.ofType(StdOutParam) == [o1]
        outs.ofType(ValueOutParam) == [v1,v2]
        outs.ofType(TupleOutParam,StdOutParam) == [o1, s1, s2]

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

    /* ------------------------------------------------
     * test path qualifier
     * ------------------------------------------------ */

    def 'test output path'() {

        given:
        def text = '''
            process foo {
              output:
              path x 
              path 'hello.*'
              path 'hello.txt'
              
              return ''
            }
            
            workflow { foo() }
            '''

        def binding = [:]

        when:
        def process = parseAndReturnProcess(text, binding)
        def out0 = (FileOutParam)process.config.getOutputs().get(0)
        def out1 = (FileOutParam)process.config.getOutputs().get(1)
        def out2 = (FileOutParam)process.config.getOutputs().get(2)

        then:
        process.config.getOutputs().size() == 3

        out0.getName() == 'x'
        out0.getFilePattern() == null
        out0.getOutChannel() instanceof DataflowVariable
        out0.isPathQualifier()

        out1.getName() == null
        out1.getFilePattern() == 'hello.*'
        out1.getOutChannel() instanceof DataflowVariable
        out1.isPathQualifier()

        out2.getName() == null
        out2.getFilePattern() == 'hello.txt'
        out2.getOutChannel() instanceof DataflowVariable
        out2.isPathQualifier()
        
    }


    def 'test output path with vars'() {

        setup:
        def text = '''

            process hola {
              output:
              path "${x}_name" 
              path "${x}_${y}.fa" 
              path "simple.txt" 
              path "data/sub/dir/file:${x}.fa" 

              return ''
            }
            
            workflow { hola() }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)
        def ctx = [x: 'hola', y:99, z:'script_file']

        when:
        FileOutParam out0 = process.config.getOutputs().get(0)
        FileOutParam out1 = process.config.getOutputs().get(1)
        FileOutParam out2 = process.config.getOutputs().get(2)
        FileOutParam out3 = process.config.getOutputs().get(3)

        then:
        process.config.getOutputs().size() == 4

        out0.name == null
        out0.getFilePatterns(ctx,null) == ['hola_name']
        out0.outChannel instanceof DataflowVariable
        out0.isDynamic()
        out0.isPathQualifier()

        out1.name == null
        out1.getFilePatterns(ctx,null) == ['hola_99.fa']
        out1.outChannel instanceof DataflowVariable
        out1.isDynamic()
        out1.isPathQualifier()

        out2.name == null
        out2.getFilePatterns(ctx,null) == ['simple.txt']
        out2.outChannel instanceof DataflowVariable
        !out2.isDynamic()
        out2.isPathQualifier()

        out3.name == null
        out3.getFilePatterns(ctx,null) == ['data/sub/dir/file:hola.fa']
        out3.outChannel instanceof DataflowVariable
        out3.isDynamic()
        out3.isPathQualifier()
    }

    def 'test output path with alias'() {

        setup:
        def text = '''

            process hola {
              output:
              path "${x}_name", emit: aaa, topic: 'foo'
              path "${x}_${y}.fa", emit: bbb 
              path "simple.txt", emit: ccc 
              path "data/sub/dir/file:${x}.fa", emit: ddd 

              return ''
            }
            
            workflow { hola() }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)
        def ctx = [x: 'hola', y:99, z:'script_file']

        when:
        FileOutParam out0 = process.config.getOutputs().get(0)
        FileOutParam out1 = process.config.getOutputs().get(1)
        FileOutParam out2 = process.config.getOutputs().get(2)
        FileOutParam out3 = process.config.getOutputs().get(3)

        then:
        process.config.getOutputs().size() == 4

        out0.name == null
        out0.getFilePatterns(ctx,null) == ['hola_name']
        out0.isDynamic()
        out0.isPathQualifier()
        out0.channelEmitName == 'aaa'
        out0.channelTopicName == 'foo'

        out1.name == null
        out1.getFilePatterns(ctx,null) == ['hola_99.fa']
        out1.outChannel instanceof DataflowVariable
        out1.isDynamic()
        out1.isPathQualifier()
        out1.channelEmitName == 'bbb'

        out2.name == null
        out2.getFilePatterns(ctx,null) == ['simple.txt']
        out2.outChannel instanceof DataflowVariable
        !out2.isDynamic()
        out2.isPathQualifier()
        out2.channelEmitName == 'ccc'

        out3.name == null
        out3.getFilePatterns(ctx,null) == ['data/sub/dir/file:hola.fa']
        out3.outChannel instanceof DataflowVariable
        out3.isDynamic()
        out3.isPathQualifier()
        out3.channelEmitName == 'ddd'
    }

    def 'test output path with tuple' () {

        given:
        def text = '''

            process hola {
              output:
              tuple path(x) 
              tuple path(y) 
              tuple path("sample.fa") 
              tuple path("data/file:${q}.fa") 

              return ''
            }
            
            workflow { hola() }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        def list_x = new BlankSeparatedList(['a.txt' as Path, 'b.txt' as Path, 'c.txt' as Path])
        def list_y = ['one.txt' as Path, 'two.txt' as Path, 'three.txt' as Path]
        def ctx = [x: list_x, y: list_y, q: 'foo']

        when:
        TupleOutParam out0 = process.config.getOutputs().get(0)
        TupleOutParam out1 = process.config.getOutputs().get(1)
        TupleOutParam out2 = process.config.getOutputs().get(2)
        TupleOutParam out3 = process.config.getOutputs().get(3)

        then:
        out0.inner[0] instanceof FileOutParam
        (out0.inner[0] as FileOutParam).getName() == 'x'
        (out0.inner[0] as FileOutParam).getFilePatterns(ctx,null) == ['a.txt', 'b.txt', 'c.txt']
        (out0.inner[0] as FileOutParam).isPathQualifier()

        out1.inner[0] instanceof FileOutParam
        (out1.inner[0] as FileOutParam).getName() == 'y'
        (out1.inner[0] as FileOutParam).getFilePatterns(ctx,null) == ['one.txt', 'two.txt', 'three.txt']
        (out1.inner[0] as FileOutParam).isPathQualifier()

        out2.inner[0] instanceof FileOutParam
        (out2.inner[0] as FileOutParam).getName() == null
        (out2.inner[0] as FileOutParam).getFilePatterns(ctx,null) == ['sample.fa']
        (out2.inner[0] as FileOutParam).isPathQualifier()

        out3.inner[0] instanceof FileOutParam
        (out3.inner[0] as FileOutParam).getName() == null
        (out3.inner[0] as FileOutParam).getFilePatterns(ctx,null) == ['data/file:foo.fa']
        (out3.inner[0] as FileOutParam).isPathQualifier()

    }

    def 'should define output path options' () {
        given:
        def text = '''

            process foo {
              output:
              path x, 
                maxDepth:2,
                hidden: false,
                followLinks: false,
                type: 'file',
                separatorChar: '#',
                glob: false,
                optional: false,
                includeInputs: false,
                arity: '1'
                
              path y, 
                maxDepth:5,
                hidden: true,
                followLinks: true,
                type: 'dir',
                separatorChar: ':',
                glob: true,
                optional: true,
                includeInputs: true,
                arity: '0..*'

              return ''
            }
            
            workflow { foo() }
            '''

        when:
        def process = parseAndReturnProcess(text, [:])
        FileOutParam out0 = process.config.getOutputs().get(0)
        FileOutParam out1 = process.config.getOutputs().get(1)

        then:
        out0.getMaxDepth() == 2
        !out0.getHidden()
        !out0.getFollowLinks()
        out0.getType()
        out0.getSeparatorChar() == '#'
        !out0.getGlob()
        !out0.getOptional()
        !out0.getIncludeInputs()
        out0.getArity() == new ArityParam.Range(1, 1)

        and:
        out1.getMaxDepth() == 5
        out1.getHidden()
        out1.getFollowLinks()
        out1.getType()
        out1.getSeparatorChar() == ':'
        out1.getGlob()
        out1.getOptional()
        out1.getIncludeInputs()
        out1.getArity() == new ArityParam.Range(0, Integer.MAX_VALUE)
    }

    def 'should set file options' () {
        given:
        def text = '''
            process foo {
              output:
              tuple path(x,maxDepth:1,optional:false), path(y,maxDepth:2,optional:true)

              return ''
            }
            
            workflow { foo() }
            '''

        when:
        def process = parseAndReturnProcess(text, [:])
        TupleOutParam out0 = process.config.getOutputs().get(0)
        then:
        out0.inner[0] instanceof FileOutParam
        and:
        (out0.inner[0] as FileOutParam).getName() == 'x'
        (out0.inner[0] as FileOutParam).getMaxDepth() == 1
        (out0.inner[0] as FileOutParam).getOptional() == false
        and:
        (out0.inner[1] as FileOutParam).getName() == 'y'
        (out0.inner[1] as FileOutParam).getMaxDepth() == 2
        (out0.inner[1] as FileOutParam).getOptional() == true
    }


    def 'should check is tuple item' () {

        setup:
        def text = '''
            
            process hola {
              output:
              val x 
              tuple val(x), path(x) 

              /command/
            }
            
            workflow { hola() }
            '''
        when:


        def process = parseAndReturnProcess(text)
        def out0 = (ValueOutParam)process.config.getOutputs().get(0)
        def out1 = (TupleOutParam)process.config.getOutputs().get(1)

        then:
        !out0.isNestedParam()
        !out1.isNestedParam()
        (out1.inner[0] as BaseOutParam).isNestedParam()
        (out1.inner[1] as BaseOutParam).isNestedParam()

        out0.index == 0
        out1.index == 1
        (out1.inner[0] as ValueOutParam).index == 1
        (out1.inner[1] as FileOutParam).index == 1
    }


    def 'should declare val with alias'() {

        setup:
        def text = '''
            process hola {
              output:
              val x,     emit: ch0
              val 10,    emit: ch1
              val 'str', emit: ch2
              val "${y}",emit: ch3
              val x.y,   emit: ch4
              /return/
            }
            
            workflow { hola() }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        def outs = process.config.getOutputs() as List<ValueOutParam>
        then:
        outs[0].name == 'x'
        outs[0].channelEmitName == 'ch0'
        outs[0].resolve([x:'foo']) == 'foo'
        and:
        outs[1].name == null
        outs[1].channelEmitName == 'ch1'
        outs[1].resolve([:]) == 10
        and:
        outs[2].name == null
        outs[2].channelEmitName == 'ch2'
        outs[2].resolve([:]) == 'str'
        and:
        outs[3].name == null
        outs[3].channelEmitName == 'ch3'
        outs[3].resolve([y:'blah']) == 'blah'
        and:
        outs[4].name == null
        outs[4].channelEmitName == 'ch4'
        outs[4].resolve( [x:[y:'Ciao']] ) == 'Ciao'

    }

    def 'should define out with emit' () {
        setup:
        def text = '''
            process hola {
              output:
              val x,     emit: ch0
              env FOO,   emit: ch1
              path '-',  emit: ch2
              stdout emit: ch3    
              /return/
            }
            
            workflow { hola() }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        def outs = process.config.getOutputs() as List<OutParam>
        then:
        outs[0].name == 'x'
        outs[0].channelEmitName == 'ch0'
        and:
        outs[1].name == 'FOO'
        outs[1].channelEmitName == 'ch1'
        and:
        outs[2] instanceof StdOutParam  // <-- note: declared as `path`, turned into a `stdout`
        outs[2].name == '-'
        outs[2].channelEmitName == 'ch2'
        and:
        outs[3] instanceof StdOutParam
        outs[3].name == '-'
        outs[3].channelEmitName == 'ch3'
    }

    def 'should define with paths' () {
        setup:
        def text = '''
            process hola {
              output:
              path x,        emit: ch0
              path 'file.*', emit: ch1
              path "${y}",   emit: ch2
              /return/
            }
            
            workflow { hola() }
            '''

        when:
        def process = parseAndReturnProcess(text)
        def outs = process.config.getOutputs() as List<FileOutParam>
        then:
        outs[0].name == 'x'
        outs[0].channelEmitName == 'ch0'
        and:
        outs[1].getName() == null
        outs[1].filePattern == 'file.*'
        outs[1].channelEmitName == 'ch1'
        and:
        outs[2].name == null
        outs[2].filePattern == null
        outs[2].channelEmitName == 'ch2'
    }

    def 'should define out tuple with alias'() {

        setup:
        def text = '''
            process hola {
              output:
                tuple val(x), val(y),   emit: ch1
                tuple path('foo'),      emit: ch2
                tuple stdout,env(bar),  emit: ch3

              /return/
            }
            
            workflow { hola() }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        def outs = process.config.getOutputs() as List<TupleOutParam>

        then:
        outs[0].name == 'tupleoutparam<0>'
        outs[0].channelEmitName == 'ch1'
        outs[0].inner[0] instanceof ValueOutParam
        outs[0].inner[0].name == 'x'
        outs[0].inner[1] instanceof ValueOutParam
        outs[0].inner[1].name == 'y'
        and:
        outs[1].name == 'tupleoutparam<1>'
        outs[1].channelEmitName == 'ch2'
        outs[1].inner[0] instanceof FileOutParam
        and:
        outs[2].name == 'tupleoutparam<2>'
        outs[2].channelEmitName == 'ch3'
        outs[2].inner[0] instanceof StdOutParam
        outs[2].inner[0].name == '-'
        outs[2].inner[1] instanceof EnvOutParam
        outs[2].inner[1].name == 'bar'

    }


    def 'should define out with topic' () {
        setup:
        def text = '''
            process hola {
              output:
              val x,     topic: ch0
              env FOO,   topic: ch1
              path '-',  topic: ch2
              stdout     topic: ch3    
              /return/
            }
            
            workflow { hola() }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        def outs = process.config.getOutputs() as List<OutParam>
        then:
        outs[0].name == 'x'
        outs[0].channelTopicName == 'ch0'
        and:
        outs[1].name == 'FOO'
        outs[1].channelTopicName == 'ch1'
        and:
        outs[2] instanceof StdOutParam  // <-- note: declared as `path`, turned into a `stdout`
        outs[2].name == '-'
        outs[2].channelTopicName == 'ch2'
        and:
        outs[3] instanceof StdOutParam
        outs[3].name == '-'
        outs[3].channelTopicName == 'ch3'
    }

    def 'should define out tuple with topic'() {

        setup:
        def text = '''
            process hola {
              output:
                tuple val(x), val(y),   topic: ch1
                tuple path('foo'),      topic: ch2
                tuple stdout,env(bar),  topic: ch3

              /return/
            }
            
            workflow { hola() }
            '''

        def binding = [:]
        def process = parseAndReturnProcess(text, binding)

        when:
        def outs = process.config.getOutputs() as List<TupleOutParam>

        then:
        outs[0].name == 'tupleoutparam<0>'
        outs[0].channelTopicName == 'ch1'
        outs[0].inner[0] instanceof ValueOutParam
        outs[0].inner[0].name == 'x'
        outs[0].inner[1] instanceof ValueOutParam
        outs[0].inner[1].name == 'y'
        and:
        outs[1].name == 'tupleoutparam<1>'
        outs[1].channelTopicName == 'ch2'
        outs[1].inner[0] instanceof FileOutParam
        and:
        outs[2].name == 'tupleoutparam<2>'
        outs[2].channelTopicName == 'ch3'
        outs[2].inner[0] instanceof StdOutParam
        outs[2].inner[0].name == '-'
        outs[2].inner[1] instanceof EnvOutParam
        outs[2].inner[1].name == 'bar'

    }
}
