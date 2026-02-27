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

package nextflow.script

import nextflow.Channel
import nextflow.exception.MissingProcessException
import nextflow.exception.ScriptCompilationException
import nextflow.exception.ScriptRuntimeException
import test.Dsl2Spec

import static test.ScriptHelper.*
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ScriptDslTest extends Dsl2Spec {

    def 'should define a process with output alias' () {
        given:
        def SCRIPT = '''

        process foo {
            output:
            val x, emit: ch1
            val y, emit: ch2

            exec:
            x = 'Hello'
            y = 'world'
        }

        workflow {
            foo()
            [ foo.out.ch1, foo.out.ch2 ]
        }
        '''

        when:
        def result = runScript(SCRIPT)
        then:
        result[0].val == 'Hello'
        result[1].val == 'world'
    }

    def 'should execute basic workflow' () {
        when:
        def result = runScript '''

        workflow hello {
            emit:
            result = 'Hello world'
        }

        workflow {
            hello()
        }
        '''

        then:
        result.val == 'Hello world'
    }

    def 'should emit expression' () {
        when:
        def result = runScript '''

        def foo() { 'Hello world' }

        workflow foo_upper {
            emit:
            foo().toUpperCase()
        }

        workflow {
            foo_upper()
        }
        '''

        then:
        result.val == 'HELLO WORLD'
    }

    def 'should emit process out' () {
        when:
        def result = runScript '''

        process foo {
          output: val x
          exec: x = 'Hello'
        }

        workflow foo_out {
            main: foo()
            emit: foo.out
        }

        workflow {
            foo_out()
        }
        '''

        then:
        result.val == 'Hello'
    }


    def 'should define processes and workflow' () {
        when:
        def result = runScript '''
        process foo {
            input: val data
            output: val result
            exec:
            result = "$data mundo"
        }

        process bar {
            input: val data
            output: val result
            exec:
            result = data.toUpperCase()
        }

        workflow alpha {
            take:
            data

            main:
            foo(data)
            bar(foo.out)

            emit:
            x = bar.out
        }

        workflow {
            alpha('Hello')
        }
        '''

        then:
        result.val == 'HELLO MUNDO'
    }


    def 'should not allow function with reserved identifier' () {

        when:
        runScript """
            def main() { println 'ciao' }
        """

        then:
        def err = thrown(ScriptCompilationException)
        err.cause.message.contains('`main` is not allowed as a function name because it is reserved for internal use')
    }

    def 'should not allow process with reserved identifier' () {

        when:
        runScript """
            process main {
              script:
              /echo ciao/
            }

            workflow {}
        """

        then:
        def err = thrown(ScriptCompilationException)
        err.cause.message.contains('`main` is not allowed as a process name because it is reserved for internal use')
    }

    def 'should not allow workflow with reserved identifier' () {

        when:
        runScript """
            workflow main {
              /echo ciao/
            }

            workflow {}
        """
        then:
        def err = thrown(ScriptCompilationException)
        err.cause.message.contains('`main` is not allowed as a workflow name because it is reserved for internal use')
    }

    def 'should not allow duplicate workflow keyword' () {
        when:
        runScript(
            """
            workflow {
              /echo ciao/
            }

            workflow {
              /echo miao/
            }
            """
        )
        then:
        def err = thrown(ScriptCompilationException)
        err.cause.message.contains('Entry workflow defined more than once')
    }

    def 'should apply operator to process result' () {
        when:
        def result = runScript('''
            process hello {
              output: val result
              exec:
                result = "Hello"
            }

            workflow {
               hello()
               hello.out.map { it.toUpperCase()  }
            }
        ''')
        then:
        result.val == 'HELLO'
    }

    def 'should branch and view' () {

        when:
        def result = runScript('''
            channel.of(1,2,3,40,50)
                .branch {
                    small: it < 10
                    large: it > 10
                }
                .set { result }

            ch1 = result.small.map { it }
            ch2 = result.large.map { it }

            [ch1, ch2]
        ''')
        then:
        result[0].val == 1
        result[0].val == 2
        result[0].val == 3
        and:
        result[1].val == 40
        result[1].val == 50
    }


    def 'should allow pipe process and operator' () {
        when:
        def result = runScript('''
        process foo {
          output: val result
          exec: result = "hello"
        }

        process bar {
          output: val result
          exec: result = "world"
        }

        workflow {
          (foo & bar) | concat
        }
        ''')

        then:
        result.val == 'hello'
        result.val == 'world'
        result.val == Channel.STOP
    }

    def 'should allow process and operator composition' () {
        when:
        def result = runScript('''
        process foo {
          output: val result
          exec: result = "hello"
        }

        process bar {
          output: val result
          exec: result = "world"
        }

        workflow {
          foo()
          bar()
          foo.out.concat(bar.out)
        }
        ''')


        then:
        result.val == 'hello'
        result.val == 'world'
        result.val == Channel.STOP
    }

    def 'should not allow invalid composition' () {
        when:
        runScript('''
        process foo {
          script:
          /echo foo/
        }

        process bar {
          input: val x
          script:
          "echo bar $x"
        }

        workflow {
          bar(foo())
        }
        ''')


        then:
        def err = thrown(ScriptRuntimeException)
        err.message == 'Process `bar` declares 1 input but was called with 0 arguments'
    }

    def 'should report error accessing undefined out/a' () {
        when:
        runScript('''
        process foo {
          script:
          /echo foo/
        }

        process bar {
          input: val x
          script:
          "echo bar $x"
        }

        workflow {
          bar(foo.out)
        }
        ''')

        then:
        def err = thrown(ScriptRuntimeException)
        err.message == "Access to 'foo.out' is undefined since the process 'foo' has not been invoked before accessing the output attribute"
    }

    def 'should report error accessing undefined out/b' () {
        when:
        runScript('''
        process foo {
          script:
          /echo foo/
        }

        process bar {
          input: val x
          script:
          "echo bar $x"
        }

        workflow {
          bar(foo.out)
        }
        ''')

        then:
        def err = thrown(ScriptRuntimeException)
        err.message == "Access to 'foo.out' is undefined since the process 'foo' has not been invoked before accessing the output attribute"
    }

    def 'should report error accessing undefined out/c' () {
        when:
        runScript('''
        process foo {
          script:
          /echo foo/
        }

        workflow flow1 {
            foo()
        }

        workflow {
          flow1()
          flow1.out.view()
        }
        ''')

        then:
        def err = thrown(ScriptRuntimeException)
        err.message == "Access to 'flow1.out' is undefined since the workflow 'flow1' doesn't declare any output"
    }

    def 'should report error accessing undefined out/d' () {
        when:
        runScript('''
        process foo {
          script:
          /echo foo/
        }

        process bar {
          input: val x
          script:
          "echo bar $x"
        }

        workflow flow1 {
            foo()
        }

        workflow {
          flow1 | bar
        }
        ''')

        then:
        def err = thrown(ScriptRuntimeException)
        err.message == "Process `bar` declares 1 input but was called with 0 arguments"
    }

    def 'should report error accessing undefined out/e' () {
        when:
        runScript('''
        process foo {
          script:
          /echo foo/
        }

        workflow flow1 {
            foo()
        }

        workflow {
          flow1.out.view()
        }
        ''')

        then:
        def err = thrown(ScriptRuntimeException)
        err.message == "Access to 'flow1.out' is undefined since the workflow 'flow1' has not been invoked before accessing the output attribute"
    }

    def 'should fail with wrong scope'() {
        when:
        runScript('''\
        process foo {
          script:
          /echo foo/
        }

        workflow foo_flow {
          main:
          flow()
          emmit:
          flow.out
        }
        ''')

        then:
        def err = thrown(ScriptCompilationException)
        err.cause.message.contains "Invalid workflow definition -- check for missing or out-of-order section labels"
    }


    def 'should fail because process is not defined'() {
        when:
        runScript(
        '''
        process sleeper {
            exec:
            """
            sleep 5
            """
        }

        workflow {
            main:
            sleeper()
            hello()
        }

        ''')

        then:
        def err = thrown(ScriptCompilationException)
        err.cause.message.contains '`hello` is not defined'
    }


    def 'should fail because is not defined /2' () {
        when:
        runScript('''
        process sleeper {
            exec:
            """
            sleep 5
            """
        }

        workflow nested {
            main:
            sleeper()
            sleeper_2()
        }

        workflow {
            nested()
        }
        ''')

        then:
        def err = thrown(ScriptCompilationException)
        err.cause.message.contains '`sleeper_2` is not defined'
    }

    def 'should fail because is not defined /3' () {
        when:
        runScript('''
        process sleeper1 {
            script:
            /echo 1/
        }

        process sleeper2 {
            script:
            /echo 3/
        }


        workflow nested {
            main:
            sleeper1()
            sleeper3()
        }

        workflow {
            nested()
        }
        ''')

        then:
        def err = thrown(ScriptCompilationException)
        err.cause.message.contains '`sleeper3` is not defined'
    }

    def 'should not conflict with private meta attribute' () {
        when:
        def result = runScript '''

        process foo {
          input: val x
          output: val y
          exec: y = x
        }

        workflow {
            meta = channel.of('Hello')
            foo(meta)
        }
        '''

        then:
        result.val == 'Hello'
    }


    def 'should throw an exception on missing method' () {

        when:
        runScript '''
           Channel.doesNotExist()
        '''
        then:
        def e1 = thrown(MissingProcessException)
        e1.message == 'Missing process or function Channel.doesNotExist()'

        when:
        runScript '''
        workflow {
           Channel.doesNotExist()
        }
        '''
        then:
        def e2 = thrown(MissingProcessException)
        e2.message == 'Missing process or function Channel.doesNotExist()'
    }

}
