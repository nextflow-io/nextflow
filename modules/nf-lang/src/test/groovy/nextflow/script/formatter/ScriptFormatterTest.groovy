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

    String format(String contents) {
        scriptParser.compiler().getSources().clear()
        def source = scriptParser.parse('main.nf', contents)
        new ScriptResolveVisitor(source, scriptParser.compiler().compilationUnit(), Types.DEFAULT_SCRIPT_IMPORTS, Collections.emptyList()).visit()
        assert !TestUtils.hasSyntaxErrors(source)
        def formatter = new ScriptFormattingVisitor(source, new FormattingOptions(4, true))
        formatter.visit()
        assert formatter.getMissingComments().isEmpty()
        return formatter.toString()
    }

    boolean checkFormat(String input, String output) {
        input = input.stripIndent()
        output = output.stripIndent()
        assert format(input) == output
        assert format(output) == output
        return true
    }

    boolean checkFormat(String source) {
        source = source.stripIndent()
        assert format(source) == source
        return true
    }

    String formatSorted(String contents) {
        scriptParser.compiler().getSources().clear()
        def source = scriptParser.parse('main.nf', contents)
        new ScriptResolveVisitor(source, scriptParser.compiler().compilationUnit(), Types.DEFAULT_SCRIPT_IMPORTS, Collections.emptyList()).visit()
        assert !TestUtils.hasSyntaxErrors(source)
        def formatter = new ScriptFormattingVisitor(source, new FormattingOptions(4, true, false, false, true))
        formatter.visit()
        assert formatter.getMissingComments().isEmpty()
        return formatter.toString()
    }

    boolean checkFormatSorted(String input, String output) {
        input = input.stripIndent()
        output = output.stripIndent()
        assert formatSorted(input) == output
        assert formatSorted(output) == output
        return true
    }

    String formatWrapped(String contents, int lineLength) {
        scriptParser.compiler().getSources().clear()
        def source = scriptParser.parse('main.nf', contents)
        new ScriptResolveVisitor(source, scriptParser.compiler().compilationUnit(), Types.DEFAULT_SCRIPT_IMPORTS, Collections.emptyList()).visit()
        assert !TestUtils.hasSyntaxErrors(source)
        def formatter = new ScriptFormattingVisitor(source, new FormattingOptions(4, true, false, false, false, lineLength))
        formatter.visit()
        assert formatter.getMissingComments().isEmpty()
        return formatter.toString()
    }

    boolean checkFormatWrapped(String input, String output, int lineLength = 40) {
        input = input.stripIndent()
        output = output.stripIndent()
        assert formatWrapped(input, lineLength) == output
        assert formatWrapped(output, lineLength) == output
        return true
    }

    boolean checkFormatWrapped(String source, int lineLength = 40) {
        source = source.stripIndent()
        assert formatWrapped(source, lineLength) == source
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

    def 'should sort includes alphabetically by source path' () {
        expect:
        checkFormatSorted(
            '''\
            include { foo } from './b.nf'
            include { bar } from './a.nf'
            ''',
            '''\
            include { bar } from './a.nf'
            include { foo } from './b.nf'
            '''
        )
    }

    def 'should leave already-sorted includes unchanged' () {
        expect:
        checkFormatSorted(
            '''\
            include { alpha } from './alpha.nf'
            include { beta } from './beta.nf'
            ''',
            '''\
            include { alpha } from './alpha.nf'
            include { beta } from './beta.nf'
            '''
        )
    }

    def 'should sort includes independently within blank-line-separated groups' () {
        expect:
        checkFormatSorted(
            '''\
            include { foo } from './b.nf'
            include { bar } from './a.nf'

            include { zed } from './z.nf'
            include { yak } from './y.nf'
            ''',
            '''\
            include { bar } from './a.nf'
            include { foo } from './b.nf'

            include { yak } from './y.nf'
            include { zed } from './z.nf'
            '''
        )
    }

    def 'should move a per-include leading comment with its include when sorted' () {
        expect:
        checkFormatSorted(
            '''\
            include { foo } from './b.nf'
            // comment about bar
            include { bar } from './a.nf'
            ''',
            '''\
            // comment about bar
            include { bar } from './a.nf'
            include { foo } from './b.nf'
            '''
        )
    }

    def 'should keep a group-header comment at the top of its group when the first include changes' () {
        expect:
        checkFormatSorted(
            '''\
            // group header
            include { foo } from './b.nf'
            include { bar } from './a.nf'
            ''',
            '''\
            // group header
            include { bar } from './a.nf'
            include { foo } from './b.nf'
            '''
        )
    }

    def 'should preserve a trailing comment on an include when sorted' () {
        expect:
        checkFormatSorted(
            '''\
            include { foo } from './b.nf' // foo comment
            include { bar } from './a.nf'
            ''',
            '''\
            include { bar } from './a.nf'
            include { foo } from './b.nf' // foo comment
            '''
        )
    }

    def 'should not sort includes when sortDeclarations is disabled' () {
        expect:
        checkFormat(
            '''\
            include { foo } from './b.nf'
            include { bar } from './a.nf'
            '''
        )
    }

    def 'should be idempotent when sorting includes across multiple groups' () {
        given:
        def output =
            '''\
            include { bar } from './a.nf'
            include { foo } from './b.nf'

            include { yak } from './y.nf'
            include { zed } from './z.nf'
            '''.stripIndent()

        expect:
        formatSorted(output) == output
        formatSorted(formatSorted(output)) == output
    }

    def 'should preserve the blank line above the include block when the first include is reordered' () {
        expect:
        checkFormatSorted(
            '''\
            nextflow.preview.output = true

            include { foo } from './b.nf'
            include { bar } from './a.nf'
            ''',
            '''\
            nextflow.preview.output = true

            include { bar } from './a.nf'
            include { foo } from './b.nf'
            '''
        )
    }

    def 'should preserve multiple blank lines between include groups' () {
        expect:
        checkFormatSorted(
            '''\
            include { foo } from './b.nf'
            include { bar } from './a.nf'


            include { zed } from './z.nf'
            include { yak } from './y.nf'
            ''',
            '''\
            include { bar } from './a.nf'
            include { foo } from './b.nf'


            include { yak } from './y.nf'
            include { zed } from './z.nf'
            '''
        )
    }

    def 'should preserve a blank line within an include leading comment block when sorted' () {
        expect:
        checkFormatSorted(
            '''\
            // header

            // note
            include { aaa } from './a.nf'
            include { bbb } from './b.nf'
            ''',
            '''\
            // header

            // note
            include { aaa } from './a.nf'
            include { bbb } from './b.nf'
            '''
        )
    }

    def 'should preserve a blank line between an include comment and its include when sorted' () {
        expect:
        checkFormatSorted(
            '''\
            include { bbb } from './b.nf'
            // note

            include { aaa } from './a.nf'
            ''',
            '''\
            // note

            include { aaa } from './a.nf'
            include { bbb } from './b.nf'
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

    // -- an `agent` declaration must survive `nextflow lint -format`. Before the AgentNode
    //    branch was added to the declaration dispatch, the node was walked and emitted NOTHING,
    //    so formatting an agent module's `main.nf` deleted its entire content.
    def 'should format an agent definition' () {
        expect:
        checkFormat(
            '''\
            nextflow.enable.types = true

            agent reporter {
                model 'openai/gpt-4o'
                instruction 'You write QA reports.'
                skills 'qa-report', 'style'

                input:
                sample: String
                depth: Integer

                output:
                report: String

                prompt:
                """
                Write a report for ${sample} at depth ${depth}.
                """
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
        checkFormat(
            '''\
            Object hello(Object x,String[] y){
            Object xy=x ; return xy
            }
            ''',
            '''\
            def hello(x, y) {
                def xy = x
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
            }
            else if (x > 0.5) {
                println('You won!')
            }
            else {
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
            }
            catch (e: IOException) {
                log.warn("Could not load foo.txt")
            }
            '''
        )
        checkFormat(
            '''\
            try{println(file('foo.txt').text)}catch(e){log.warn("Could not load foo.txt")}
            ''',
            '''\
            try {
                println(file('foo.txt').text)
            }
            catch (e) {
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

    /// COMMENTS

    def 'should preserve a trailing comment' () {
        expect:
        checkFormat(
            '''\
            params.index = null     // path to the index
            ''',
            '''\
            params.index = null // path to the index
            '''
        )
        checkFormat(
            '''\
            workflow {
                x = 1 // about x
                y = 2 // about y
            }
            '''
        )
        checkFormat(
            '''\
            workflow {
                x = 1 // the last statement
            }
            '''
        )
    }

    def 'should preserve a dangling comment' () {
        expect:
        // at the end of a block
        checkFormat(
            '''\
            workflow {
                x = channel.of(1)

                // PROCESS_1(x)
                // PROCESS_2(x)
            }
            '''
        )
        // in an empty block
        checkFormat(
            '''\
            workflow {
                // nothing here yet
            }
            '''
        )
        // at the end of the file
        checkFormat(
            '''\
            workflow {
                println('hi')
            }

            // process something {
            // }
            '''
        )
    }

    def 'should preserve comments in a list or map literal' () {
        expect:
        checkFormat(
            '''\
            workflow {
                x = [
                    // first
                    1,
                    2, // second
                ]
                y = [
                    a: 1, // one
                    // two
                    b: 2,
                ]
            }
            '''
        )
    }

    def 'should preserve comments in call arguments' () {
        expect:
        checkFormat(
            '''\
            workflow {
                FOO(
                    // first arg
                    1,
                    2, // second arg
                )
            }
            '''
        )
    }

    def 'should preserve comments in a process' () {
        expect:
        checkFormat(
            '''\
            process FOO {
                // a directive
                cpus 2

                input:
                // about x
                val x

                output:
                path 'out.txt' // the output

                script:
                """
                echo hi
                """
            }
            '''
        )
    }

    def 'should preserve comments in an agent' () {
        expect:
        checkFormat(
            '''\
            nextflow.enable.types = true

            agent reporter {
                // a directive
                model 'openai/gpt-4o'

                input:
                // about sample
                sample: String

                output:
                report: String // the report

                prompt:
                """
                Write a report for ${sample}.
                """
            }
            '''
        )
    }

    def 'should preserve comments in an if/else' () {
        expect:
        checkFormat(
            '''\
            workflow {
                if (true) {
                    // then branch
                    x = 1
                }
                else {
                    // else branch
                    x = 2
                }
            }
            '''
        )
    }

    def 'should preserve a comment before an else branch' () {
        expect:
        checkFormat(
            '''\
            workflow {
                if (x) {
                    a()
                }
                // otherwise
                else if (y) {
                    b()
                }
                // last resort
                else {
                    c()
                }
            }
            '''
        )
    }

    def 'should preserve a block comment' () {
        expect:
        checkFormat(
            '''\
            /*
             * A block comment
             */
            workflow {
                println('hi')
            }
            '''
        )
    }

    def 'should preserve a file header' () {
        expect:
        checkFormat(
            '''\
            #!/usr/bin/env nextflow
            /*
             * Copyright 2013-2024, Seqera Labs
             */

            // enable the output DSL
            nextflow.preview.output = true

            /*
             * Command line input parameter
             */
            params.query = 'data/sample.fa'
            '''
        )
    }

    def 'should preserve comments in a params block' () {
        expect:
        checkFormat(
            '''\
            params {
                // Comma-separated list of IDs.
                input: String

                // Whether to save intermediate outputs.
                save_intermeds: Boolean = false
            }

            workflow {
                println(params.input)
            }
            '''
        )
    }

    def 'should preserve comments on a function' () {
        expect:
        checkFormat(
            '''\
            /*
             * does a thing
             */
            def foo(x) {
                return x // give it back
            }
            '''
        )
    }

    def 'should preserve a trailing comment on a single-line if' () {
        expect:
        // the branch is expanded onto its own line, so the comment trails
        // the statement as a whole
        checkFormat(
            '''\
            workflow {
                if (!group) group = [1] // if new, create a new entry
            }
            ''',
            '''\
            workflow {
                if (!group) {
                    group = [1]
                } // if new, create a new entry
            }
            '''
        )
    }

    def 'should preserve a comment after a closure parameter list' () {
        expect:
        checkFormat(
            '''\
            workflow {
                x = [1, 2, 2, 3].reduce([:]) { acc, v -> // 'acc' collects all values
                    acc[v] = 1
                    return acc // give it back
                }
            }
            ''',
            '''\
            workflow {
                x = [1, 2, 2, 3].reduce([:]) { acc, v ->
                    // 'acc' collects all values
                    acc[v] = 1
                    return acc // give it back
                }
            }
            '''
        )
    }

    def 'should preserve comments in a crlf source' () {
        given:
        def source = "// lead\r\nworkflow {\r\n    x = 1   // trail\r\n}\r\n"

        expect:
        format(source) == "// lead\nworkflow {\n    x = 1 // trail\n}\n"
    }

    def 'should hoist a comment that has no place to go' () {
        expect:
        checkFormat(
            '''\
            workflow {
                x = 1 +
                    // about the second operand
                    2
            }
            ''',
            '''\
            workflow {
                // about the second operand
                x = 1 + 2
            }
            '''
        )
    }

    def 'should preserve comments on method chain links' () {
        expect:
        // an end-of-line comment trails the link it follows
        checkFormat(
            '''\
            workflow {
                result = channel.of(1, 2)
                    .map { it * 2 } // double
                    .filter { it > 2 }
                    .view()
            }
            '''
        )
        // an own-line comment leads the link it precedes
        checkFormat(
            '''\
            workflow {
                result = channel.of(1, 2)
                    .map { it * 2 }
                    // about the filter
                    .filter { it > 2 }
                    .view()
            }
            '''
        )
        // a comment forces the chain to wrap even when the source is one line
        checkFormat(
            '''\
            workflow {
                result = channel.of(1, 2).map { it * 2 } // double
                    .view()
            }
            ''',
            '''\
            workflow {
                result = channel.of(1, 2)
                    .map { it * 2 } // double
                    .view()
            }
            '''
        )
        // a namespace receiver is broken onto its own line for a comment
        checkFormat(
            '''\
            workflow {
                channel // about of
                    .of(1)
                    .view()
            }
            '''
        )
    }

    /// MULTI-LINE STRINGS

    def 'should re-indent a multi-line string with its statement' () {
        expect:
        // the statement moves right, so the string body moves right
        checkFormat(
            '''\
            process FOO {
              script:
              """
              echo hi
                echo there
              """
            }
            ''',
            '''\
            process FOO {
                script:
                """
                echo hi
                  echo there
                """
            }
            '''
        )
        // the statement moves left, so the string body moves left
        checkFormat(
            '''\
            process FOO {
                script:
                  """
                  echo hi
                    echo there
                  """
            }
            ''',
            '''\
            process FOO {
                script:
                """
                echo hi
                  echo there
                """
            }
            '''
        )
    }

    def 'should never shift a string line past its content' () {
        expect:
        // `beta` has no indentation to give up, so it stays put and the
        // relative indentation of the string body is not preserved
        checkFormat(
            '''\
            workflow {
                    x = """
                      alpha
            beta
                      """
                    println(x)
            }
            ''',
            '''\
            workflow {
                x = """
                  alpha
            beta
                  """
                println(x)
            }
            '''
        )
    }

    def 'should not re-indent a string whose statement starts mid-line' () {
        expect:
        // the closure body moved because the chain was wrapped, not because
        // its indentation changed, so there is nothing meaningful to shift by
        checkFormat(
            '''\
            workflow {
              channel.of(1).map { v -> v }.view { """
              value: ${v}
              """ }
            }
            ''',
            '''\
            workflow {
                channel.of(1)
                    .map { v -> v }
                    .view {
                        """
              value: ${v}
              """
                    }
            }
            '''
        )
    }

    // -- line-length wrapping (issue #26)

    def 'should wrap a long function call and leave a short one alone' () {
        expect:
        checkFormatWrapped(
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
            ''',
            60
        )
        checkFormatWrapped(
            '''\
            workflow {
                x = foo(a, b)
            }
            ''',
            '''\
            workflow {
                x = foo(a, b)
            }
            ''',
            60
        )
    }

    def 'should wrap a long list literal and leave a short one alone' () {
        expect:
        checkFormatWrapped(
            '''\
            workflow {
                x = [alpha, beta, gamma, delta, epsilon, zeta]
            }
            ''',
            '''\
            workflow {
                x = [
                    alpha,
                    beta,
                    gamma,
                    delta,
                    epsilon,
                    zeta,
                ]
            }
            '''
        )
        checkFormatWrapped(
            '''\
            workflow {
                x = [alpha, beta, gamma]
            }
            ''',
            '''\
            workflow {
                x = [alpha, beta, gamma]
            }
            '''
        )
    }

    def 'should wrap a long map literal' () {
        expect:
        checkFormatWrapped(
            '''\
            workflow {
                x = [alpha: 1, beta: 2, gamma: 3, delta: 4]
            }
            ''',
            '''\
            workflow {
                x = [
                    alpha: 1,
                    beta: 2,
                    gamma: 3,
                    delta: 4,
                ]
            }
            '''
        )
    }

    def 'should wrap a long process directive' () {
        expect:
        checkFormatWrapped(
            '''\
            process FOO {
                publishDir params.outdir, mode: 'copy', overwrite: true

                script:
                """
                echo hi
                """
            }
            ''',
            '''\
            process FOO {
                publishDir(
                    params.outdir,
                    mode: 'copy',
                    overwrite: true,
                )

                script:
                """
                echo hi
                """
            }
            '''
        )
    }

    def 'should wrap a long structured process input' () {
        expect:
        checkFormatWrapped(
            '''\
            nextflow.enable.types = true

            process FOO {
                input:
                tuple(identifier: String, someInputFile: Path)

                script:
                'echo hi'
            }
            ''',
            '''\
            nextflow.enable.types = true

            process FOO {
                input:
                tuple(
                    identifier: String,
                    someInputFile: Path
                )

                script:
                'echo hi'
            }
            '''
        )
    }

    def 'should escalate to wrapping every construct when the outer wrap alone is not enough' () {
        expect:
        checkFormatWrapped(
            '''\
            workflow {
                x = foo(alpha, [nested_one, nested_two, nested_three, nested_four])
            }
            ''',
            '''\
            workflow {
                x = foo(
                    alpha,
                    [
                        nested_one,
                        nested_two,
                        nested_three,
                        nested_four,
                    ],
                )
            }
            '''
        )
    }

    def 'should preserve a comment inside a force-wrapped list' () {
        expect:
        checkFormatWrapped(
            '''\
            workflow {
                x = [alpha, beta, gamma, delta, epsilon, zeta] // the letters
            }
            ''',
            '''\
            workflow {
                x = [
                    alpha,
                    beta,
                    gamma,
                    delta,
                    epsilon,
                    zeta,
                ] // the letters
            }
            '''
        )
    }

    def 'should not wrap inside a string interpolation' () {
        expect:
        checkFormatWrapped(
            '''\
            process FOO {
                script:
                "echo ${[alpha, beta, gamma].join(' ')} and some more text past the limit"
            }
            '''
        )
    }

    def 'should not wrap a long line when maxLineLength is disabled' () {
        expect:
        checkFormatWrapped(
            '''\
            workflow {
                ALIGN_AND_SORT(samples_channel, reference_genome, annotation_file, params.threads)
            }
            ''',
            '''\
            workflow {
                ALIGN_AND_SORT(samples_channel, reference_genome, annotation_file, params.threads)
            }
            ''',
            0
        )
    }

}
