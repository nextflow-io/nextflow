/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.script.formatter

import nextflow.script.control.ScriptParser
import nextflow.script.control.ScriptResolveVisitor
import nextflow.script.dsl.Types
import spock.lang.Shared
import spock.lang.Specification
import test.TestUtils

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptFormatterTest extends Specification {

    @Shared
    ScriptParser scriptParser

    def setupSpec() {
        scriptParser = new ScriptParser()
    }

    String format(FormattingOptions options, String contents) {
        scriptParser.compiler().getSources().clear()
        def source = scriptParser.parse('main.nf', contents)
        new ScriptResolveVisitor(source, scriptParser.compiler().compilationUnit(), Types.DEFAULT_SCRIPT_IMPORTS, Collections.emptyList()).visit()
        assert !TestUtils.hasSyntaxErrors(source)
        def formatter = new ScriptFormattingVisitor(source, options, contents)
        formatter.visit()
        return formatter.toString()
    }

    String format(String contents) {
        return format(new FormattingOptions(4, true), contents)
    }

    boolean checkFormat(FormattingOptions options, String input, String output) {
        input = input.stripIndent()
        output = output.stripIndent()
        assert format(options, input) == output
        assert format(options, output) == output
        return true
    }

    boolean checkFormat(String input, String output) {
        return checkFormat(new FormattingOptions(4, true), input, output)
    }

    boolean checkFormat(String source) {
        source = source.stripIndent()
        assert format(source) == source
        return true
    }

    boolean checkFormat(FormattingOptions options, String source) {
        source = source.stripIndent()
        assert format(options, source) == source
        return true
    }

    def 'should format a code snippet' () {
        expect:
        checkFormat(
            '''\
            println 'Hello!'
            ''',
            '''\
            println('Hello!')
            '''
        )
    }

    /// SCRIPT DECLARATIONS

    def 'should format an include declaration' () {
        expect:
        checkFormat(
            '''\
            include{foo;bar}from'./foobar.nf'
            ''',
            '''\
            include { foo ; bar } from './foobar.nf'
            '''
        )
        checkFormat(
            '''\
            include{
            foo;bar
            }from'./foobar.nf'
            ''',
            '''\
            include {
                foo ;
                bar
            } from './foobar.nf'
            '''
        )
    }

    def 'should format a parameter declaration' () {
        expect:
        checkFormat(
            '''\
            params.foo='bar'
            ''',
            '''\
            params.foo = 'bar'
            '''
        )
    }

    def 'should format a legacy workflow definition' () {
        expect:
        checkFormat(
            '''\
            workflow hello{
            take: x ; y ; main: xy=x*y ; emit: result = xy
            }
            ''',
            '''\
            workflow hello {
                take:
                x
                y

                main:
                xy = x * y

                emit:
                result = xy
            }
            '''
        )

        checkFormat(
            '''\
            workflow hello{
            take: x ; y ; emit: result = x * y
            }
            ''',
            '''\
            workflow hello {
                take:
                x
                y

                emit:
                result = x * y
            }
            '''
        )
    }

    def 'should format a typed workflow definition' () {
        expect:
        checkFormat(
            '''\
            nextflow.enable.types = true

            workflow hello{
            take: x:Integer ; y:Integer ; main: xy=x*y ; emit: result:Integer = xy
            }
            ''',
            '''\
            nextflow.enable.types = true

            workflow hello {
                take:
                x: Integer
                y: Integer

                main:
                xy = x * y

                emit:
                result: Integer = xy
            }
            '''
        )
    }

    def 'should format a legacy process definition' () {
        expect:
        checkFormat(
            '''\
            process hello{
            debug(true) ; input: val x ; val y ; output: tuple val(x), val(y) ; script: 'echo true'
            }
            ''',
            '''\
            process hello {
                debug true

                input:
                val x
                val y

                output:
                tuple val(x), val(y)

                script:
                'echo true'
            }
            '''
        )
    }

    def 'should format a typed process definition' () {
        expect:
        checkFormat(
            '''\
            nextflow.enable.types=true

            process hello{
            debug(true) ; input: tuple(id:String,infile:Path) ; index:Path ; stage: stageAs(infile,'input.txt') ; output: result=tuple(id,file('output.txt')) ; script: 'cat input.txt > output.txt'
            }
            ''',
            '''\
            nextflow.enable.types = true

            process hello {
                debug true

                input:
                tuple(id: String, infile: Path)
                index: Path

                stage:
                stageAs infile, 'input.txt'

                output:
                result = tuple(id, file('output.txt'))

                script:
                'cat input.txt > output.txt'
            }
            '''
        )

        checkFormat(
            '''\
            nextflow.enable.types=true

            process hello{
            input: record(id:String,infile:Path) ; script: 'cat input.txt > output.txt'
            }
            ''',
            '''\
            nextflow.enable.types = true

            process hello {
                input:
                record(
                    id: String,
                    infile: Path
                )

                script:
                'cat input.txt > output.txt'
            }
            '''
        )
    }

    def 'should format a function definition' () {
        expect:
        checkFormat(
            '''\
            def hello(x,y){
            def xy=x*y ; return xy
            }
            ''',
            '''\
            def hello(x, y) {
                def xy = x * y
                return xy
            }
            '''
        )
        checkFormat(
            '''\
            Integer hello(Integer x,Integer y){
            Integer xy=x*y ; return xy
            }
            ''',
            '''\
            def hello(x: Integer, y: Integer) -> Integer {
                def xy: Integer = x * y
                return xy
            }
            '''
        )
        checkFormat(
            '''\
            def hello(x:Integer,y:Integer)->Integer{
            def xy=x*y ; return xy
            }
            ''',
            '''\
            def hello(x: Integer, y: Integer) -> Integer {
                def xy = x * y
                return xy
            }
            '''
        )
    }

    def 'should format an enum definition' () {
        expect:
        checkFormat(
            '''\
            enum Colors{RED,GREEN,BLUE}
            ''',
            '''\
            enum Colors {
                RED,
                GREEN,
                BLUE,
            }
            '''
        )
    }

    def 'should format a record definition' () {
        expect:
        checkFormat(
            '''\
            record FastqPair{id:String;fastq_1: Path;fastq_2: Path?}
            ''',
            '''\
            record FastqPair {
                id: String
                fastq_1: Path
                fastq_2: Path?
            }
            '''
        )
    }

    def 'should format an output block' () {
        expect:
        checkFormat(
            '''\
            workflow{}

            output{
            foo{path'foo'}
            bar{path'bar';index{path'index.csv'}}
            }
            ''',
            '''\
            workflow {
            }

            output {
                foo {
                    path 'foo'
                }
                bar {
                    path 'bar'
                    index {
                        path 'index.csv'
                    }
                }
            }
            '''
        )
        checkFormat(
            '''\
            workflow{}

            output{
            foo:Path{path'foo'}
            bar:Channel<Path>{path'bar';index{path'index.csv'}}
            }
            ''',
            '''\
            workflow {
            }

            output {
                foo: Path {
                    path 'foo'
                }
                bar: Channel<Path> {
                    path 'bar'
                    index {
                        path 'index.csv'
                    }
                }
            }
            '''
        )
    }

    def 'should not sort script declarations by default' () {
        given:
        def source =
            '''\
            params.foo = 'bar'

            include { foo ; bar } from './foobar.nf'

            process hello {
                script:
                'echo true'
            }

            workflow hello {
            }

            workflow {
            }

            output {
            }
            '''.stripIndent()

        expect:
        format(source) == source
    }

    /// STATEMENTS

    def 'should format an assert statement' () {
        expect:
        checkFormat(
            '''\
            assert 2+2==4:'The math broke!'
            ''',
            '''\
            assert 2 + 2 == 4 : 'The math broke!'
            '''
        )
    }

    def 'should format a variable declaration' () {
        expect:
        checkFormat(
            '''\
            def x=42
            def(x,y)=tuple(1,2)
            def(
            x,
            y
            )=tuple(1,2)
            ''',
            '''\
            def x = 42
            def (x, y) = tuple(1, 2)
            def (x, y) = tuple(1, 2)
            '''
        )
        checkFormat(
            '''\
            def Integer x=42
            ''',
            '''\
            def x: Integer = 42
            '''
        )
        checkFormat(
            '''\
            def x:Integer=42
            def(x:Integer,y:Integer)=tuple(1,2)
            ''',
            '''\
            def x: Integer = 42
            def (x: Integer, y: Integer) = tuple(1, 2)
            '''
        )
    }

    def 'should format an assignment' () {
        expect:
        checkFormat(
            '''\
            v=42
            list[0]='first'
            map.key='value'
            (x,y)=tuple(1,2)
            (
            x,
            y
            )=tuple(1,2)
            ''',
            '''\
            v = 42
            list[0] = 'first'
            map.key = 'value'
            (x, y) = tuple(1, 2)
            (x, y) = tuple(1, 2)
            '''
        )
    }

    def 'should format an if-else statement' () {
        expect:
        checkFormat(
            '''\
            if(x<0.5)println('You lost.')else if(x>0.5)println('You won!')else println('You tied?')
            ''',
            '''\
            if (x < 0.5) {
                println('You lost.')
            } else if (x > 0.5) {
                println('You won!')
            } else {
                println('You tied?')
            }
            '''
        )
    }

    def 'should format a throw statement' () {
        expect:
        checkFormat(
            '''\
            throw new Exception('something failed!')
            '''
        )
    }

    def 'should format a try-catch statement' () {
        expect:
        checkFormat(
            '''\
            try{println(file('foo.txt').text)}catch(IOException e){log.warn("Could not load foo.txt")}
            ''',
            '''\
            try {
                println(file('foo.txt').text)
            } catch (e: IOException) {
                log.warn("Could not load foo.txt")
            }
            '''
        )
    }

    def 'should format a method chain' () {
        expect:
        checkFormat(
            '''\
            channel.of( 1, 2, 3 )
                .multiMap{v->foo:bar:v}.set{result}
            ''',
            '''\
            channel.of(1, 2, 3)
                .multiMap { v -> foo: bar: v }
                .set { result }
            '''
        )
    }

    /// EXPRESSIONS

    def 'should format a numeric literal' () {
        expect:
        checkFormat(
            '''\
            42
            -1
            0b1001
            031
            0xabcd
            3.14
            -0.1
            1.59e7
            1.59e-7
            '''
        )
    }

    def 'should format a keyword literal' () {
        expect:
        checkFormat(
            '''\
            true
            false
            null
            '''
        )
    }

    def 'should format a string literal' () {
        expect:
        checkFormat(
            '''\
            "I said 'hello'"

            'I said "hello" again!'

            \'\'\'
            Hello,
            How are you today?
            \'\'\'

            """
            We don't have to escape quotes anymore!
            Even "double" quotes!
            """

            /no escape!/
            '''
        )
    }

    def 'should format a dynamic string' () {
        expect:
        checkFormat(
            '''\
            "Hello, ${names.join(' and ')}!"
            '''
        )
        checkFormat(
            '''\
            "Hello, $name.first $name.last!"
            ''',
            '''\
            "Hello, ${name.first} ${name.last}!"
            '''
        )
        checkFormat(
            '''\
            """
            blastp \
                -in ${input} \
                -out ${output} \
                -db ${blast_db} \
                -html
            """
            '''
        )
        checkFormat(
            '''\
            'Hello, ${names.join(" and ")}!'
            '''
        )
    }

    def 'should format a list literal' () {
        expect:
        checkFormat(
            '''\
            [1,2,3]
            []
            ''',
            '''\
            [1, 2, 3]
            []
            '''
        )
    }

    def 'should format a map literal' () {
        expect:
        checkFormat(
            '''\
            [foo:1,bar:2,baz:3]
            [:]
            ''',
            '''\
            [foo: 1, bar: 2, baz: 3]
            [:]
            '''
        )
        checkFormat(
            '''\
            [(x):1]
            ''',
            '''\
            [(x): 1]
            '''
        )
        checkFormat(
            '''\
            [(x.y):1]
            ''',
            '''\
            [(x.y): 1]
            '''
        )
    }

    def 'should format a closure' () {
        expect:
        checkFormat(
            '''\
            {a,b->a+b}
            ''',
            '''\
            { a, b -> a + b }
            '''
        )
        checkFormat(
            '''\
            {Integer a,Integer b->a+b}
            ''',
            '''\
            { a: Integer, b: Integer -> a + b }
            '''
        )
        checkFormat(
            '''\
            {a:Integer,b:Integer->a+b}
            ''',
            '''\
            { a: Integer, b: Integer -> a + b }
            '''
        )
        checkFormat(
            '''\
            {v->println'Hello!';v*v}
            ''',
            '''\
            { v ->
                println('Hello!')
                v * v
            }
            '''
        )
    }

    def 'should format index and property accesses' () {
        expect:
        checkFormat(
            '''\
            myList[0]
            myFile.text
            myFiles*.text
            myFile?.text
            '''
        )
    }

    def 'should format a function call' () {
        expect:
        checkFormat(
            '''\
            printf 'Hello %s!\\n', 'World'
            file 'hello.txt', checkIfExists: true
            ''',
            '''\
            printf('Hello %s!\\n', 'World')
            file('hello.txt', checkIfExists: true)
            '''
        )
        checkFormat(
            '''\
            [1, 2, 3].inject('result:') { acc, v -> acc + ' ' + v }
            [1, 2, 3].each() { v -> println v }
            [1, 2, 3].each { v -> println v }
            ''',
            '''\
            [1, 2, 3].inject('result:') { acc, v -> acc + ' ' + v }
            [1, 2, 3].each { v -> println(v) }
            [1, 2, 3].each { v -> println(v) }
            '''
        )
    }

    def 'should format a constructor call' () {
        expect:
        checkFormat(
            '''\
            new java.util.Date()
            new Date()
            '''
        )
    }

    def 'should format unary/binary/ternary expressions' () {
        expect:
        checkFormat(
            '''\
            !(2+2==4)
            (1+2)*3
            x%2==0?'x is even!':'x is odd!'
            (false?'foo':true)?'bar':'baz'
            ''',
            '''\
            !(2 + 2 == 4)
            (1 + 2) * 3
            x % 2 == 0 ? 'x is even!' : 'x is odd!'
            (false ? 'foo' : true) ? 'bar' : 'baz'
            '''
        )
    }

    // -- comment preservation (issues #111, #127, #140)

    def 'should preserve comments at the end of the file' () {
        expect:
        // https://github.com/nextflow-io/language-server/issues/111
        checkFormat(
            '''\
            workflow {
                Channel.of("a", "b", "c", "d") | view
            }

            // process something {
            //   input:
            //   val a
            //
            //   script:
            //   """
            //   touch b.txt
            //   """
            // }
            ''',
            '''\
            workflow {
                Channel.of("a", "b", "c", "d") | view
            }

            // process something {
            //   input:
            //   val a
            //
            //   script:
            //   """
            //   touch b.txt
            //   """
            // }
            '''
        )
    }

    def 'should preserve comments at the end of a block' () {
        expect:
        // https://github.com/nextflow-io/language-server/issues/140
        checkFormat(
            '''\
            workflow {
                fastq_files_ch.view()

                // PROCESS_1(fastq_files_ch)
                // new_ch = PROCESS_1.out
                // PROCESS_2(new_ch)
            }
            ''',
            '''\
            workflow {
                fastq_files_ch.view()

                // PROCESS_1(fastq_files_ch)
                // new_ch = PROCESS_1.out
                // PROCESS_2(new_ch)
            }
            '''
        )
        checkFormat(
            '''\
            def foo() {
                return 42
                // trailing note
            }

            workflow {
                if( params.x ) {
                    println 'hi'
                    // done greeting
                }
                else {
                    // nothing to do
                }
            }
            ''',
            '''\
            def foo() {
                return 42
                // trailing note
            }

            workflow {
                if (params.x) {
                    println('hi')
                    // done greeting
                } else {
                    // nothing to do
                }
            }
            '''
        )
    }

    def 'should preserve comments in an empty block' () {
        expect:
        checkFormat(
            '''\
            workflow {
                // TODO: implement
            }
            ''',
            '''\
            workflow {
                // TODO: implement
            }
            '''
        )
    }

    def 'should keep trailing comments on the same line' () {
        expect:
        // https://github.com/nextflow-io/language-server/issues/140
        checkFormat(
            '''\
            params.index_ref = null     // Path to index reference genome
            params.other = 1

            workflow {
                x = 1 // set x
                y = 2 // set y
            }
            ''',
            '''\
            params.index_ref = null // Path to index reference genome
            params.other = 1

            workflow {
                x = 1 // set x
                y = 2 // set y
            }
            '''
        )
    }

    def 'should format leading comments in a file with CRLF line endings' () {
        expect:
        // https://github.com/nextflow-io/language-server/issues/127
        checkFormat(
            "workflow {\r\n    // ALIGN reads to reference genome\r\n    BWA_ALIGN(sample_id, library_id)\r\n}\r\n",
            '''\
            workflow {
                // ALIGN reads to reference genome
                BWA_ALIGN(sample_id, library_id)
            }
            '''
        )
    }

    def 'should preserve comments in process sections' () {
        expect:
        checkFormat(
            '''\
            process foo {
                tag "sample"
                // this used to be cpus 4

                input:
                // the sample id
                val sample_id

                output:
                // the result
                path "result.txt" // inline comment

                script:
                """
                touch result.txt
                """
                // after script
            }
            ''',
            '''\
            process foo {
                tag "sample"
                // this used to be cpus 4

                input:
                // the sample id
                val sample_id

                output:
                // the result
                path "result.txt" // inline comment

                script:
                """
                touch result.txt
                """
                // after script
            }
            '''
        )
    }

    def 'should preserve comments in workflow sections' () {
        expect:
        checkFormat(
            '''\
            workflow FOO {
                take:
                // the samples
                samples // inline take
                // after takes

                main:
                result = samples

                emit:
                // the first output
                out1 = result
                // the second output
                out2 = result // inline emit
                // after emits
            }
            ''',
            '''\
            workflow FOO {
                take:
                // the samples
                samples // inline take
                // after takes

                main:
                result = samples

                emit:
                // the first output
                out1 = result
                // the second output
                out2 = result // inline emit
                // after emits
            }
            '''
        )
    }

    def 'should preserve comments in a params block' () {
        expect:
        checkFormat(
            '''\
            params {
                // the input file
                input: Path
                // an option
                opt: String = 'default' // inline note
                // dangling params comment
            }

            workflow {
                println params.input
            }
            ''',
            '''\
            params {
                // the input file
                input: Path
                // an option
                opt: String = 'default' // inline note
                // dangling params comment
            }

            workflow {
                println(params.input)
            }
            '''
        )
    }

    def 'should preserve comments in closures' () {
        expect:
        checkFormat(
            '''\
            workflow {
                ch = Channel.of(1, 2).map { v ->
                    v + 1
                    // add one
                }
            }
            ''',
            '''\
            workflow {
                ch = Channel
                    .of(1, 2)
                    .map { v ->
                        v + 1
                        // add one
                    }
            }
            '''
        )
    }

    def 'should keep comments in place inside wrapped expressions' () {
        expect:
        // comments before the links of a method chain stay at their links
        checkFormat(
            '''\
            workflow {
                ch = Channel.of(1, 2)
                    // filter odd numbers
                    .filter { v -> v % 2 == 0 }
                    // double the values
                    .map { v -> v * 2 }
            }
            ''',
            '''\
            workflow {
                ch = Channel
                    .of(1, 2)
                    // filter odd numbers
                    .filter { v -> v % 2 == 0 }
                    // double the values
                    .map { v -> v * 2 }
            }
            '''
        )
        // comments before the elements of a wrapped call or collection stay
        // at their elements
        checkFormat(
            '''\
            workflow {
                PINTS_CALLER(
                    // the grouped bams
                    group_bam_bai,
                    // the assay
                    params.assay_type,
                )
                chromosomes = [
                    // autosomes only
                    "chr1",
                    "chr2"]
            }
            ''',
            '''\
            workflow {
                PINTS_CALLER(
                    // the grouped bams
                    group_bam_bai,
                    // the assay
                    params.assay_type,
                )
                chromosomes = [
                    // autosomes only
                    "chr1",
                    "chr2",
                ]
            }
            '''
        )
        // comments in positions the formatter cannot emit in place are
        // hoisted above the statement instead of being removed
        checkFormat(
            '''\
            workflow {
                x = (1 +
                    // half of it
                    2)
            }
            ''',
            '''\
            workflow {
                // half of it
                x = (1 + 2)
            }
            '''
        )
    }

    // -- K&R style (issue #153)

    def 'should format if-else and try-catch in K&R style' () {
        expect:
        // https://github.com/nextflow-io/language-server/issues/153
        checkFormat(
            '''\
            workflow {
                if( params.a ) {
                    run_a()
                }
                else if( params.b ) {
                    run_b()
                }
                else {
                    run_c()
                }
            }
            ''',
            '''\
            workflow {
                if (params.a) {
                    run_a()
                } else if (params.b) {
                    run_b()
                } else {
                    run_c()
                }
            }
            '''
        )
        checkFormat(
            '''\
            def f() {
                try {
                    g()
                }
                catch( Exception e ) {
                    h()
                }
            }

            workflow {
                f()
            }
            ''',
            '''\
            def f() {
                try {
                    g()
                } catch (e: Exception) {
                    h()
                }
            }

            workflow {
                f()
            }
            '''
        )
        // a comment above a catch clause stays above it (the clause falls
        // back to its own line so the comment can be emitted)
        checkFormat(
            '''\
            def f() {
                try {
                    g()
                }
                // when things go wrong
                catch (e: Exception) {
                    h()
                }
            }

            workflow {
                f()
            }
            ''',
            '''\
            def f() {
                try {
                    g()
                }
                // when things go wrong
                catch (e: Exception) {
                    h()
                }
            }

            workflow {
                f()
            }
            '''
        )
    }

    // -- blank line normalization (issues #115, #150)

    def 'should normalize blank lines' () {
        expect:
        // https://github.com/nextflow-io/language-server/issues/115
        // https://github.com/nextflow-io/language-server/issues/150
        checkFormat(
            '''\
            include { A } from './a.nf'
            include { B } from './b.nf'
            params.x = 1
            workflow FOO {
                foo()
            }
            process bar {
                script:
                "true"
            }
            workflow {
                FOO()
            }
            ''',
            '''\
            include { A } from './a.nf'
            include { B } from './b.nf'

            params.x = 1

            workflow FOO {
                foo()
            }

            process bar {
                script:
                "true"
            }

            workflow {
                FOO()
            }
            '''
        )
        checkFormat(
            '''\
            workflow {

                x = 1



                y = 2

            }
            ''',
            '''\
            workflow {
                x = 1

                y = 2
            }
            '''
        )
        checkFormat(
            '''\
            #!/usr/bin/env nextflow
            workflow {
                x = 1
            }
            ''',
            '''\
            #!/usr/bin/env nextflow

            workflow {
                x = 1
            }
            '''
        )
    }

    // -- multi-line string re-indentation (issue #116)

    def 'should re-indent multi-line strings' () {
        expect:
        // https://github.com/nextflow-io/language-server/issues/116
        checkFormat(
            new FormattingOptions(2, true, false, false, false, 120),
            '''\
            process foo {
                script:
                """
                echo 'hello world!'
                """
            }
            ''',
            '''\
            process foo {
              script:
              """
              echo 'hello world!'
              """
            }
            '''
        )
        checkFormat(
            '''\
            workflow {
            x = """
                line one
                  line two
                """
            }
            ''',
            '''\
            workflow {
                x = """
                    line one
                      line two
                    """
            }
            '''
        )
    }

    // -- fmt: skip / fmt: off / fmt: on (issue #75)

    def 'should not format regions excluded with fmt directives' () {
        expect:
        // https://github.com/nextflow-io/language-server/issues/75
        checkFormat(
            '''\
            workflow {
                x  =  [1,  2,   3] // fmt: skip
                y = [4,5]
            }
            ''',
            '''\
            workflow {
                x  =  [1,  2,   3] // fmt: skip
                y = [4, 5]
            }
            '''
        )
        checkFormat(
            '''\
            workflow {
                a = 1

                // fmt: off
                matrix = [
                    [1, 0],
                    [0, 1] ]
                // fmt: on

                b = 2
            }
            ''',
            '''\
            workflow {
                a = 1

                // fmt: off
                matrix = [
                    [1, 0],
                    [0, 1] ]
                // fmt: on

                b = 2
            }
            '''
        )
    }

    // -- line-length wrapping (issue #26)

    def 'should wrap lines that exceed the maximum line length' () {
        given:
        def options = new FormattingOptions(4, true, false, false, false, 60)

        expect:
        // https://github.com/nextflow-io/language-server/issues/26
        checkFormat(options,
            '''\
            workflow {
                ALIGN_AND_SORT(samples_channel, reference_genome, annotation_file, params.threads)
            }
            ''',
            '''\
            workflow {
                ALIGN_AND_SORT(
                    samples_channel,
                    reference_genome,
                    annotation_file,
                    params.threads,
                )
            }
            '''
        )
        checkFormat(options,
            '''\
            workflow {
                result = Channel.fromPath(params.input).splitCsv(header: true).map { row -> row.sample }
            }
            ''',
            '''\
            workflow {
                result = Channel
                    .fromPath(params.input)
                    .splitCsv(header: true)
                    .map { row -> row.sample }
            }
            '''
        )
        // short lines are not wrapped
        checkFormat(options,
            '''\
            workflow {
                x = foo(a, b)
            }
            ''',
            '''\
            workflow {
                x = foo(a, b)
            }
            '''
        )
    }

    // -- include sorting (issue #54)

    def 'should sort includes within groups when sorting is enabled' () {
        given:
        def options = new FormattingOptions(4, true, false, false, true, 120)

        expect:
        // https://github.com/nextflow-io/language-server/issues/54
        checkFormat(options,
            '''\
            // local modules
            include { ZULU } from './modules/zulu.nf'
            include { ALPHA } from './modules/alpha.nf'

            // subworkflows
            include { UTILS } from '../subworkflows/utils.nf'
            include { BAM_SORT } from '../subworkflows/bam_sort.nf'

            workflow {
                ALPHA()
            }
            ''',
            '''\
            // local modules
            include { ALPHA } from './modules/alpha.nf'
            include { ZULU } from './modules/zulu.nf'

            // subworkflows
            include { BAM_SORT } from '../subworkflows/bam_sort.nf'
            include { UTILS } from '../subworkflows/utils.nf'

            workflow {
                ALPHA()
            }
            '''
        )
    }

    def 'should keep trailing comments on wrapped elements and chain links' () {
        expect:
        // a comment on the same line as an element or chain link stays on
        // that line
        checkFormat(
            '''\
            workflow {
                ALIGN(
                    reads, // the raw reads
                    genome, // the reference
                )
                ch = Channel.of(1, 2) // make a channel
                    .map { v -> v * 2 } // double it
                    .view()
            }
            ''',
            '''\
            workflow {
                ALIGN(
                    reads, // the raw reads
                    genome, // the reference
                )
                ch = Channel
                    .of(1, 2) // make a channel
                    .map { v -> v * 2 } // double it
                    .view()
            }
            '''
        )
    }

    def 'should keep comments before the closing bracket of a construct' () {
        expect:
        checkFormat(
            '''\
            workflow {
                samples = [
                    "a",
                    "b",
                    // more to come
                ]
                ALIGN(
                    reads,
                    genome,
                    // add extra args here
                )
            }
            ''',
            '''\
            workflow {
                samples = [
                    "a",
                    "b",
                    // more to come
                ]
                ALIGN(
                    reads,
                    genome,
                    // add extra args here
                )
            }
            '''
        )
    }

    def 'should keep comments on the branches of a wrapped ternary' () {
        expect:
        checkFormat(
            '''\
            workflow {
                x = params.aligner == 'bwa' // check the aligner
                    // use the fast path
                    ? 'fast'
                    // otherwise be safe
                    : 'safe'
            }
            ''',
            '''\
            workflow {
                x = params.aligner == 'bwa' // check the aligner
                    // use the fast path
                    ? 'fast'
                    // otherwise be safe
                    : 'safe'
            }
            '''
        )
    }

    def 'should keep blank-line grouping around element comments' () {
        expect:
        checkFormat(
            '''\
            workflow {
                samples = [
                    "a",
                    "b",

                    // the odd one
                    "z",
                ]
            }
            ''',
            '''\
            workflow {
                samples = [
                    "a",
                    "b",

                    // the odd one
                    "z",
                ]
            }
            '''
        )
    }

    def 'should keep chain comments in section and param values' () {
        expect:
        // method chains in emit values, params and config values are
        // wrapped so their comments stay at their links
        checkFormat(
            '''\
            params {
                input: String = 'default'
                    // strip whitespace
                    .trim()
            }

            workflow W {
                main:
                x = Channel.of(1)

                emit:
                out = x
                    // squared values
                    .map { v -> v * v }
            }

            workflow {
                W()
            }
            ''',
            '''\
            params {
                input: String = 'default'
                    // strip whitespace
                    .trim()
            }

            workflow W {
                main:
                x = Channel.of(1)

                emit:
                out = x
                    // squared values
                    .map { v -> v * v }
            }

            workflow {
                W()
            }
            '''
        )
        checkFormat(
            '''\
            params.outdir = file(params.base)
                // resolve the results dir
                .resolve('results')

            workflow {
                println params.outdir
            }
            ''',
            '''\
            params.outdir = file(params.base)
                // resolve the results dir
                .resolve('results')

            workflow {
                println(params.outdir)
            }
            '''
        )
    }

    def 'should preserve comments in a code snippet' () {
        expect:
        checkFormat(
            '''\
            // emit some letters
            channel.of('a'..'z') | view

            // emit some numbers
            channel.of(1..9) | view
            ''',
            '''\
            // emit some letters
            channel.of('a'..'z') | view

            // emit some numbers
            channel.of(1..9) | view
            '''
        )
    }

    def 'should preserve blank lines after compound statements' () {
        expect:
        checkFormat(
            '''\
            def rule(file) {
                if (file == 'file_1.txt') {
                    return "alpha/${file}"
                }

                if (file == 'file_2.txt') {
                    return null
                }
            }

            workflow {
                rule('file_1.txt')
            }
            '''
        )
    }

    def 'should keep a method chain with property links wrapped' () {
        expect:
        checkFormat(
            '''\
            workflow {
                foo()
                foo.out
                    .collect()
                    .view { v -> "got: ${v}" }
            }
            '''
        )
    }

    // -- workflow completion handlers

    def 'should emit the main label when a workflow has completion handlers' () {
        expect:
        checkFormat(
            '''\
            workflow {
                main:
                foo()

                onComplete:
                println('done')
            }
            '''
        )
        checkFormat(
            '''\
            workflow {
                main:
                foo()

                onError:
                println('failed')
            }
            '''
        )
    }

    // -- trailing comments after compound statements

    def 'should keep a trailing comment after a compound statement' () {
        expect:
        checkFormat(
            '''\
            workflow {
                if (params.check) {
                    foo()
                } // done checking
                bar()
            }
            '''
        )
        checkFormat(
            '''\
            def f() {
                try {
                    g()
                } catch (e: Exception) {
                    h()
                } // recovered
                i()
            }

            workflow {
                f()
            }
            '''
        )
    }

    // -- fmt directives on section labels, directives and fields

    def 'should not format directives and fields excluded with fmt directives' () {
        expect:
        // without the directive, the process directive is formatted
        checkFormat(
            '''\
            process foo {
                cpus 4*2

                script:
                'true'
            }
            ''',
            '''\
            process foo {
                cpus 4 * 2

                script:
                'true'
            }
            '''
        )
        // process directive
        checkFormat(
            '''\
            process foo {
                cpus 4*2 // fmt: skip

                input:
                val x

                script:
                'true'
            }
            '''
        )
        // record field
        checkFormat(
            '''\
            nextflow.enable.types = true

            record Point {
                x : Integer // fmt: skip
                y: Integer
            }

            workflow {
                println('ok')
            }
            '''
        )
        // when section
        checkFormat(
            '''\
            process foo {
                input:
                val x

                when:
                x    <    1 // fmt: skip

                script:
                'true'
            }
            '''
        )
    }

    // -- string escapes (issues #41, #55)

    def 'should preserve string escapes' () {
        expect:
        checkFormat(
            '''\
            workflow {
                a = "quote \\" and backslash \\\\ and dollar \\$HOME"
                b = 'single \\' quote'
                c = "keep ${x} interp"
            }
            '''
        )
    }

    // -- spread operator (issue #70)

    def 'should preserve the spread operator in method calls' () {
        expect:
        checkFormat(
            '''\
            workflow {
                ch.map { it -> tuple(it[0], it[1]*.get(0)) }.view()
            }
            '''
        )
    }

    // -- magic trailing comma (issues #60, #74)

    def 'should keep source-wrapped calls wrapped and preserve trailing commas' () {
        expect:
        checkFormat(
            '''\
            workflow {
                foo(
                    a,
                    b
                )
                bar(
                    a,
                    b,
                )
            }
            ''',
            '''\
            workflow {
                foo(
                    a,
                    b,
                )
                bar(
                    a,
                    b,
                )
            }
            '''
        )
    }

    // -- line wrapping can be disabled

    def 'should not wrap lines when the maximum line length is disabled' () {
        expect:
        checkFormat(
            new FormattingOptions(4, true, false, false, false, 0),
            '''\
            workflow {
                ALIGN_AND_SORT(samples_channel, reference_genome, annotation_file, params.threads, params.aligner_options)
            }
            ''',
            '''\
            workflow {
                ALIGN_AND_SORT(samples_channel, reference_genome, annotation_file, params.threads, params.aligner_options)
            }
            '''
        )
    }

    // -- harshil alignment (issues #71, #135)

    def 'should align includes and emits with harshil alignment' () {
        given:
        def options = new FormattingOptions(4, true, true, false, false, 120)

        expect:
        checkFormat(options,
            '''\
            include { FOO } from './foo.nf'
            include { LONGER_NAME } from './modules.nf'
            ''',
            '''\
            include { FOO         } from './foo.nf'
            include { LONGER_NAME } from './modules.nf'
            '''
        )
        checkFormat(options,
            '''\
            workflow FOO {
                main:
                x = 1
                result = 2

                emit:
                out = x
                longer_name = result
            }
            ''',
            '''\
            workflow FOO {
                main:
                x = 1
                result = 2

                emit:
                out         = x
                longer_name = result
            }
            '''
        )
        // legacy process outputs are not aligned (see language-server #71):
        // typed processes don't have this problem
        checkFormat(options,
            '''\
            process bar {
                input:
                val x

                output:
                tuple val(x), path("*.txt"), emit: results
                path "versions.yml", emit: versions

                script:
                'true'
            }
            '''
        )
    }

    // -- mahesh form (issue #129)

    def 'should move process outputs after the script with mahesh form' () {
        expect:
        checkFormat(
            new FormattingOptions(4, true, false, true, false, 120),
            '''\
            process foo {
                input:
                val x

                output:
                path 'y.txt'

                script:
                'true'
            }
            ''',
            '''\
            process foo {
                input:
                val x

                script:
                'true'

                output:
                path 'y.txt'
            }
            '''
        )
    }

    // -- comment safety check API

    def 'should collect comment texts for the safety check' () {
        expect:
        CommentReattacher.countComments('// a\nx = 1 /* b */\n', false) == 2
        CommentReattacher.commentTexts('// a\nx = 1 /* b */\n', false) == ['/* b */', '// a']
        // CRLF line endings are normalized
        CommentReattacher.commentTexts('// a\r\nx = 1 /* b\r\n c */\r\n', false) == ['/* b\n c */', '// a']
        // the shebang is not counted as a comment
        CommentReattacher.countComments('#!/usr/bin/env nextflow\n// a\n', false) == 1
        CommentReattacher.commentTexts('#!/usr/bin/env nextflow\n// a\n', false) == ['// a']
    }

}
