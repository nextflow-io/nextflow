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

import java.nio.file.Files

import nextflow.exception.MissingProcessException
import nextflow.exception.ScriptCompilationException
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.*
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class ScriptIncludesTest extends Dsl2Spec {

    def 'should invoke foreign functions' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        def alpha() {
          return 'this is alpha result'
        }

        def bravo(x) {
          return x.reverse()
        }

        def gamma(x,y) {
          return "$x and $y"
        }
        '''

        SCRIPT.text = """
        include { alpha; bravo; gamma } from "$MODULE"

        def local_func() {
          return "I'm local"
        }

        workflow {
          [
            alpha(),
            bravo('Hello'),
            gamma('Hola', 'mundo'),
            local_func()
          ]
        }
        """

        when:
        def result = runScript(SCRIPT)

        then:
        result[0] == 'this is alpha result'
        result[1] == 'olleH'
        result[2] == 'Hola and mundo'
        result[3] == "I'm local"
    }

    def 'should invoke foreign functions from operator' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        def foo(str) {
          return str.reverse()
        }
        '''

        SCRIPT.text = """
        include { foo } from "$MODULE"
        workflow {
          channel.value('hello world').map { foo(it) }
        }
        """

        when:
        def result = runScript(SCRIPT)

        then:
        result.val == 'dlrow olleh'
    }

    def 'should allow functions with default args' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        def foo(str='foo') {
          return str.reverse()
        }
        '''

        SCRIPT.text = """
        include { foo } from "$MODULE"
        workflow {
          [ witharg: foo('hello world'), withdefault: foo() ]
        }
        """

        when:
        def result = runScript(SCRIPT)
        then:
        result instanceof Map
        result.witharg == 'hello world'.reverse()
        result.withdefault == 'foo'.reverse()
    }

    def 'should allow multiple signatures of function' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        def foo( list=[1,2,3] ) {
          return list
        }
        def foo(c1, c2){
            return c1+"-"+c2
        }
        '''

        SCRIPT.text = """
        include { foo } from "$MODULE"
        workflow {
          foo().collect { foo(it, it*2) }
        }
        """

        when:
        def result = runScript(SCRIPT)

        then:
        result[0] == '1-2'
        result[1] == '2-4'
        result[2] == '3-6'
    }

    def 'should fail if no signatures of function founded' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        def foo( list=[1,2,3] ) {
          return list
        }
        def foo(c1, c2){
            return c1+"-"+c2
        }
        '''

        SCRIPT.text = """
        include { foo } from "$MODULE"
        workflow {
          channel.of( foo(1, 2, 3) )
        }
        """

        when:
        runScript(SCRIPT)

        then:
        thrown(MissingProcessException)
    }

    def 'should invoke a workflow from include' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
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
            take: data
            main: foo(data)
                  bar(foo.output)
            emit: bar.out
        }
        '''

        SCRIPT.text = """
        include { alpha } from "$MODULE"

        workflow {
            alpha('Hello')
        }
        """

        when:
        def result = runScript(SCRIPT)

        then:
        result.val == 'HELLO MUNDO'
    }


    def 'should invoke a workflow from main' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
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
        '''

        SCRIPT.text = """
        include { foo; bar } from "$MODULE"

        workflow alpha {
            take: data
            main: foo(data)
                  bar(foo.output)
            emit: bar.out
        }

        workflow {
            alpha('Hello')
        }
        """

        when:
        def result = runScript(SCRIPT)

        then:
        result.val == 'HELLO MUNDO'

    }

    def 'should invoke an entry workflow' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
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
        '''

        SCRIPT.text = """
        include { foo; bar } from "$MODULE"

        workflow {
            data = 'Hello'
            bar(foo(data))
        }
        """

        when:
        def result = runScript(SCRIPT)

        then:
        result.val == 'HELLO MUNDO'
    }

    def 'should gather process outputs' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
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
        '''

        SCRIPT.text = """
        include { foo; bar } from "$MODULE"

        workflow {
            data = 'Hello'
            foo(data)
            bar(foo.output)
            bar.out
        }
        """

        when:
        def result = runScript(SCRIPT)

        then:
        result.val == 'HELLO MUNDO'
    }

    def 'should define a process and invoke it' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        process foo {
          input: val sample
          output: stdout
          script:
          "echo Hello $sample"
        }
        '''

        SCRIPT.text = """
        include { foo } from "$MODULE"

        workflow {
          foo('world')
        }
        """

        when:
        def result = runScript(SCRIPT)
        then:
        noExceptionThrown()
        result.val == 'echo Hello world'

        cleanup:
        folder?.deleteDir()
    }

    def 'should define a process with multiple inputs' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        process foo {
          input:
            val sample
            tuple val(pairId), val(reads)
          output:
            stdout
          script:
            "echo sample=$sample pairId=$pairId reads=$reads"
        }
        '''

        SCRIPT.text = """
        include { foo } from './module.nf'

        workflow {
          ch1 = channel.of('world')
          ch2 = channel.value(['x', '/some/file'])
          foo(ch1, ch2)
        }
        """

        when:
        def result = runScript(SCRIPT)
        then:
        noExceptionThrown()
        result.val == 'echo sample=world pairId=x reads=/some/file'

        cleanup:
        folder?.deleteDir()
    }


    def 'should compose processes' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        process foo {
          input:
            val alpha
          output:
            val delta
            val gamma
          script:
            delta = alpha
            gamma = 'world'
            /nope/
        }

        process bar {
           input:
             val xx
             val yy
           output:
             stdout
           script:
             "echo $xx $yy"
        }
        '''

        and:
        SCRIPT.text = """
        include { foo; bar } from './module.nf'

        workflow {
            (delta, gamma) = foo('Ciao')
            bar( delta, gamma )
        }
        """

        when:
        def result = runScript(SCRIPT)
        then:
        result.val == 'echo Ciao world'

        cleanup:
        folder?.deleteDir()
    }

    def 'should use multiple assignment' () {

        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        process foo {
          input:
            val alpha
          output:
            val delta
            val gamma
          script:
            delta = alpha
            gamma = 'world'
            /nope/
        }

        process bar {
           input:
             val xx
             val yy
           output:
             stdout
           script:
             "echo $xx $yy"
        }
        '''

        and:
        SCRIPT.text = """
        include { foo } from './module.nf'

        workflow {
          (ch0, ch1) = foo('Ciao')
          [ ch0, ch1 ]
        }
        """

        when:
        def result = runScript(SCRIPT)
        then:
        result[0].val == 'Ciao'
        result[1].val == 'world'

        cleanup:
        folder?.deleteDir()
    }


    def 'should invoke custom functions' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        def foo(str) {
          str.reverse()
        }

        def bar(a, b) {
          return "$a $b!"
        }
        '''

        SCRIPT.text = """
        include { foo; bar } from './module.nf'

        workflow {
          def str = foo('dlrow')
          return bar('Hello', str)
        }
        """

        when:
        def result = runScript(SCRIPT)
        then:
        noExceptionThrown()
        result == 'Hello world!'

        cleanup:
        folder?.deleteDir()
    }

    def 'should not fail when invoking a process in a module' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        process foo {
            script:
            /hello/
        }

        workflow { foo() }
        '''

        SCRIPT.text = """
        include { foo } from './module.nf'

        workflow {
            println 'hello'
        }
        """

        when:
        runScript(SCRIPT)
        then:
        noExceptionThrown()

        cleanup:
        folder?.deleteDir()
    }

    def 'should include modules' () {
        given:
        def folder = Files.createTempDirectory('test')
        folder.resolve('org').mkdir()
        def MODULE = folder.resolve('org/bio.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        process foo {
            script:
            /hello/
        }
        '''

        SCRIPT.text = """
        include { foo } from './org/bio'

        workflow {
            foo()
        }
        """

        when:
        runScript(SCRIPT)
        then:
        noExceptionThrown()

        cleanup:
        folder?.deleteDir()
    }

    def 'should include only named component' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        def alpha() {
          return 'this is alpha result'
        }

        def bravo(x) {
          return x.reverse()
        }
        '''

        when:
        SCRIPT.text = """
        include { alpha } from "$MODULE"

        workflow {
          alpha()
        }
        """
        def result = runScript(SCRIPT)
        then:
        result == 'this is alpha result'

        when:
        SCRIPT.text = """
        include { alpha as FOO } from "$MODULE"

        workflow {
          FOO()
        }
        """
        result = runScript(SCRIPT)
        then:
        result == 'this is alpha result'


        when:
        SCRIPT.text = """
        include { alpha as FOO } from "$MODULE"

        workflow {
          alpha()
        }
        """
        result = runScript(SCRIPT)
        then:
        thrown(ScriptCompilationException)


        when:
        SCRIPT.text = """
        include { alpha } from "$MODULE"

        workflow {
          bravo()
        }
        """
        result = runScript(SCRIPT)
        then:
        thrown(ScriptCompilationException)

        cleanup:
        folder?.deleteDir()
    }


    def 'should invoke process aliases' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        process foo {
          input: val data
          output: val result
          exec:
            result = data.toUpperCase()
        }
        '''

        SCRIPT.text = """
        include { foo } from "$MODULE"
        include { foo as bar } from "$MODULE"

        workflow {
            [ foo('Hello'), bar('World') ]
        }
        """

        when:
        def result = runScript(SCRIPT)
        then:
        result[0].val == 'HELLO'
        result[1].val == 'WORLD'

        cleanup:
        folder?.deleteDir()
    }


    def 'should allow multiple use of the same process' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
            process producer {
                output: stdout
                script: "echo Hello"
            }

            process consumer {
                input: file "foo"
                output: stdout
                script:
                "cmd consumer 1"
            }

            process another_consumer {
                input: file "foo"
                output: stdout
                script: "cmd consumer 2"
            }

            workflow flow1 {
                emit: producer | consumer | map { it.toUpperCase() }
            }

            workflow flow2 {
                emit: producer | another_consumer | map { it.toUpperCase() }
            }
            '''.stripIndent()

        SCRIPT.text = """
            include { flow1; flow2 } from "$MODULE"

            workflow {
              flow1()
              flow2()

              [ flow1.out, flow2.out ]
            }

            """

        when:
        def result = runScript(SCRIPT)

        then:
        result[0].val == 'CMD CONSUMER 1'
        result[1].val == 'CMD CONSUMER 2'

        cleanup:
        folder?.deleteDir()
    }



    def 'should include many processes' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        process foo {
          input: val data
          output: val result
          exec:
            result = data.toUpperCase()
        }
        '''

        SCRIPT.text = """
        include { foo; foo as bar } from "$MODULE"

        workflow {
            foo('Hello')
            bar('World')
            [ foo.out, bar.out ]
        }
        """

        when:
        def result = runScript(SCRIPT)

        then:
        result[0].val == 'HELLO'
        result[1].val == 'WORLD'

        cleanup:
        folder?.deleteDir()
    }

    def 'should declare moduleDir path' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module/dir').createDirIfNotExists().resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = """
        def foo () { return true }

        workflow {
          // moduleDir uses real path (symlinks resolved) for consistent lookups
          assert moduleDir == file("$folder/module/dir")
          assert projectDir == file("$folder")
          assert launchDir == file('.')
        }
        """

        SCRIPT.text = """
        include { foo } from "$MODULE"

        workflow {
          // moduleDir uses real path (symlinks resolved) for consistent lookups
          assert moduleDir == file("$folder")
          assert projectDir == file("$folder")
          assert launchDir == file('.')
        }
        """

        when:
        runScript(SCRIPT)
        then:
        noExceptionThrown()

        cleanup:
        folder?.deleteDir()
    }

    def 'should not allow include nested within a workflow' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''

        process foo {
          script:
          /echo hello/
        }
        '''

        SCRIPT.text = """
        workflow {
            include { foo } from "./module.nf"
            foo()
        }
        """

        when:
        runScript(SCRIPT)
        then:
        thrown(ScriptCompilationException)

        cleanup:
        folder?.deleteDir()
    }

    def 'should should allow invoking function passing gstring' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        def alpha(String str) {
            return str.reverse()
        }
        '''

        when:
        SCRIPT.text = """
        include { alpha } from "$MODULE"

        workflow {
            def x = "world"
            def y = "Hello \$x"

            alpha(y)
        }
        """
        def result = runScript(SCRIPT)
        then:
        result == 'dlrow olleH'

        cleanup:
        folder?.deleteDir()
    }
}
