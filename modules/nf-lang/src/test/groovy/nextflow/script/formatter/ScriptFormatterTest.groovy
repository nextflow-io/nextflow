/*
 * Copyright 2024-2025, Seqera Labs
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
import nextflow.script.types.Types
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

    def 'should format a workflow definition' () {
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

        checkFormat(
            '''\
            workflow hello{
            take: x:Integer ; y:Integer ; main: xy=x*y ; emit: result:Integer = xy
            }
            ''',
            '''\
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
            nextflow.preview.types=true

            process hello{
            debug(true) ; input: (id,infile):Tuple<String,Path> ; index:Path ; stage: stageAs('input.txt',infile) ; output: result=tuple(id,file('output.txt')) ; script: 'cat input.txt > output.txt'
            }
            ''',
            '''\
            nextflow.preview.types = true

            process hello {
                debug true

                input:
                (id, infile): Tuple<String, Path>
                index: Path

                stage:
                stageAs 'input.txt', infile

                output:
                result = tuple(id, file('output.txt'))

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

    def 'should format an output block' () {
        expect:
        checkFormat(
            '''\
            output{
            foo{path'foo'}
            bar{path'bar';index{path'index.csv'}}
            }
            ''',
            '''\
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
            output{
            foo:Path{path'foo'}
            bar:Channel<Path>{path'bar';index{path'index.csv'}}
            }
            ''',
            '''\
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
            ''',
            '''\
            !(2 + 2 == 4)
            (1 + 2) * 3
            x % 2 == 0 ? 'x is even!' : 'x is odd!'
            '''
        )
    }

}
