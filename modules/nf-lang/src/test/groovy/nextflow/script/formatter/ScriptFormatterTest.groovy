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
        def source = scriptParser.parse('main.nf', contents)
        new ScriptResolveVisitor(source, scriptParser.compiler().compilationUnit(), Types.DEFAULT_IMPORTS, Collections.emptyList()).visit()
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
    }

    def 'should format a process definition' () {
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
    }

    def 'should format an operator chain' () {
        expect:
        checkFormat(
            '''\
            Channel.of( 1, 2, 3 )
                .multiMap{v->foo:bar:v}.set{result}
            ''',
            '''\
            Channel.of(1, 2, 3)
                .multiMap { v -> foo: bar: v }
                .set { result }
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

}
