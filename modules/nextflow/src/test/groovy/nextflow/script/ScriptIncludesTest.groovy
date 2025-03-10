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

import nextflow.NextflowMeta
import nextflow.exception.MissingProcessException
import nextflow.exception.ScriptCompilationException
import spock.lang.Timeout
import test.Dsl2Spec
import test.MockScriptRunner
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class ScriptIncludesTest extends Dsl2Spec {

    def 'should catch wrong script' () {
        given:
        def test = Files.createTempDirectory('test')
        def lib = Files.createDirectory(test.toAbsolutePath()+"/lib")
        def MODULE = lib.resolve('Foo.groovy')
        def SCRIPT = test.resolve('main.nf')

        MODULE.text = '''
        class Foo {
            String id
        }
        '''

        SCRIPT.text = """
        include { Foo } from "$MODULE" 
        
        process foo {
            input:
                val value
        
            output:
                path '*.txt'
        
            script:
                "echo 'hello'"
        }
        workflow {
            foo(Channel.of(new Foo(id: "hello_world")))
        }
        """

        when:
        new MockScriptRunner().setScript(SCRIPT).execute()

        then:
        def err = thrown(ScriptCompilationException)
        err.message == """\
                Module compilation error
                - file : $MODULE
                """.stripIndent().rightTrim()
    }

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
   
        ret1 = alpha()
        ret2 = bravo('Hello')
        ret3 = gamma('Hola', 'mundo')
        ret4 = local_func()
        """

        when:
        def runner = new MockScriptRunner()
        def binding = runner.session.binding
        def result = runner.setScript(SCRIPT).execute()

        then:
        binding.ret1 == 'this is alpha result'
        binding.ret2 == 'olleH'
        binding.ret3 == 'Hola and mundo'
        binding.ret4 == "I'm local"
        binding.ret4 == result
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
           emit:
           channel.of('hello world').map { foo(it) }
        }
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result.val == 'dlrow olleh'
    }

    def 'should allow duplicate functions' () {
        given:
        NextflowMeta.instance.strictMode(true)
        and:
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
           emit:
           channel.of('hello world').map { 
            [ witharg : foo(it), withdefault : foo() ] 
           }
        }
        """

        when:
        def result = new MockScriptRunner() .setScript(SCRIPT).execute()
        def map = result.val
        then:
        map
        map.witharg == 'hello world'.reverse()
        map.withdefault == 'foo'.reverse()
        
        cleanup:
        NextflowMeta.instance.strictMode(false)
    }

    def 'should allow multiple signatures of function' () {
        given:
        NextflowMeta.instance.strictMode(true)
        and:
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
           emit:
           channel.fromList( foo() ).flatMap { foo(it, it*2) } 
        }
        """

        when:
        def result = new MockScriptRunner() .setScript(SCRIPT).execute()

        then:
        result.val == '1-2'
        result.val == '2-4'
        result.val == '3-6'

        cleanup:
        NextflowMeta.instance.strictMode(false)
    }

    def 'should fail if no signatures of function founded' () {
        given:
        NextflowMeta.instance.strictMode(true)
        and:
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
           emit:
           channel.of( foo(1, 2, 3) ) 
        }
        """

        when:
        def result = new MockScriptRunner() .setScript(SCRIPT).execute()

        then:
        thrown(MissingProcessException)

        cleanup:
        NextflowMeta.instance.strictMode(false)
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
            main: alpha('Hello')
            emit: alpha.out 
        }
        """

        when:
        def runner = new MockScriptRunner()
        def binding = runner.session.binding
        def result = runner.setScript(SCRIPT).execute()

        then:
        result.val == 'HELLO MUNDO'
        binding.variables.alpha == null
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
            main: alpha('Hello')
            emit: alpha.out 
        }
        """

        when:
        def runner = new MockScriptRunner()
        def binding = runner.session.binding
        def result = runner.setScript(SCRIPT).execute()

        then:
        result.val == 'HELLO MUNDO'
        !binding.hasVariable('alpha')
        !binding.hasVariable('foo')
        !binding.hasVariable('bar')

    }

    def 'should invoke an anonymous workflow' () {
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
   
        data = 'Hello'
        workflow {
            main: foo(data)
                  bar(foo.output)
            emit: bar.out 
        }
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

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
            emit: bar.out 
        }
        """

        when:
        def runner = new MockScriptRunner()
        def vars = runner.session.binding.variables
        def result = runner.setScript(SCRIPT).execute()

        then:
        result.val == 'HELLO MUNDO'
        !vars.containsKey('data')
        !vars.containsKey('foo')
        !vars.containsKey('bar')
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
          /echo Hello $sample/
        }        
        '''

        SCRIPT.text = """
        include { foo } from "$MODULE" 
        hello_ch = Channel.of('world')
        
        workflow {
            main: foo(hello_ch)
            emit: foo.out 
        }
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        noExceptionThrown()
        result.val == 'echo Hello world'

        cleanup:
        folder?.deleteDir()
    }

    def 'should define a process with multiple inputs' () {
        given:
        def folder = TestHelper.createInMemTempDir()
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
            /echo sample=$sample pairId=$pairId reads=$reads/
        }
        '''

        SCRIPT.text = """
        include { foo } from './module.nf'

        workflow {
          main: ch1 = Channel.of('world')
                ch2 = Channel.value(['x', '/some/file'])
                foo(ch1, ch2)
          emit: foo.out  
        }
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        noExceptionThrown()
        result.val == 'echo sample=world pairId=x reads=/some/file'
    }


    def 'should compose processes' () {
        given:
        def folder = TestHelper.createInMemTempDir()
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
            /echo $xx $yy/            
        }
        '''

        and:
        SCRIPT.text = """  
        include { foo; bar } from './module.nf'        

        workflow {
            main: bar( foo('Ciao') )
            emit: bar.out
        }
        """

        when:
        def result = dsl_eval(SCRIPT)
        then:
        result.val == 'echo Ciao world'
    }

    def 'should use multiple assignment' () {

        given:
        def folder = TestHelper.createInMemTempDir()
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
            /echo $xx $yy/            
        }
        '''

        and:
        SCRIPT.text = """ 
        include { foo } from './module.nf'        
        
        workflow {
          main: (ch0, ch1) = foo('Ciao')
          emit: ch0; ch1
        }
        """
        
        when:
        def result = dsl_eval(SCRIPT)
        then:
        result[0].val == 'Ciao'
        result[1].val == 'world'
    }


    def 'should inject params in module' () {
        given:
        def folder = TestHelper.createInMemTempDir()
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        params.foo = 'x' 
        params.bar = 'y'
        
        process foo {
          output: stdout 
          script:
          /echo $params.foo $params.bar/
        }
        '''

        // inject params in the module
        // and invoke the process 'foo'
        SCRIPT.text = """     
        include { foo } from "./module.nf" params(foo:'Hello', bar: 'world')
            
        workflow { 
            main: foo()
            emit: foo.out 
        }
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        noExceptionThrown()
        result.val == 'echo Hello world'
        
    }


    def 'should invoke custom functions' () {
        given:
        def folder = TestHelper.createInMemTempDir()
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

        def str = foo('dlrow')
        return bar('Hello', str)
        """

        when:
        def result = new MockScriptRunner()
                .setScript(SCRIPT)
                .execute()
        then:
        noExceptionThrown()
        result == 'Hello world!'
    }

    def 'should access module variables' () {
        given:
        def folder = TestHelper.createInMemTempDir()
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''     
        params.x = 'Hello world'
        FOO = params.x   
        
        process foo {
          output: stdout 
          script:
          "echo $FOO"
        }
        '''

        SCRIPT.text = """ 
        include { foo } from './module.nf' params(x: 'Hola mundo')
        
        workflow {
            main: foo()
            emit: foo.out
        }    
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        noExceptionThrown()
        result.val == 'echo Hola mundo'
    }

    def 'should not fail when invoking a process in a module' () {
        given:
        def folder = TestHelper.createInMemTempDir()
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''     
        process foo {
            /hello/ 
        }      
        
        workflow { foo() } 
        '''

        SCRIPT.text = """ 
        include { foo } from './module.nf'
        println 'hello'
        """

        when:
        def runner = new MockScriptRunner()
        runner.setScript(SCRIPT).execute()
        then:
        noExceptionThrown()
    }

    def 'should include modules' () {
        given:
        def folder = TestHelper.createInMemTempDir();
        folder.resolve('org').mkdir()
        def MODULE = folder.resolve('org/bio.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''     
        process foo {
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
        def runner = new MockScriptRunner()
        runner.setScript(SCRIPT).execute()
        then:
        noExceptionThrown()
    }

    def 'should include only named component' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')

        MODULE.text = '''
        def alpha() {
          return 'this is alpha result'
        }   
        
        def bravo(x) { 
          return x.reverse()
        }

        '''

        when:
        def SCRIPT = """
        include { alpha } from "$MODULE"
        alpha()
        """
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        result == 'this is alpha result'

        when:
        SCRIPT = """
        include { alpha as FOO } from "$MODULE"
        FOO()
        """
        runner = new MockScriptRunner()
        result = runner.setScript(SCRIPT).execute()
        then:
        result == 'this is alpha result'


        when:
        SCRIPT = """
        include { alpha as FOO } from "$MODULE"
        alpha()
        """
        runner = new MockScriptRunner()
        runner.setScript(SCRIPT).execute()
        then:
        thrown(MissingMethodException)


        when:
        SCRIPT = """
        include { alpha } from "$MODULE"
        bravo()
        """
        runner = new MockScriptRunner()
        runner.setScript(SCRIPT).execute()
        then:
        thrown(MissingMethodException)

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
            foo('Hello')
            bar('World')
            emit: foo.out
            emit: bar.out
        }
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result[0].val == 'HELLO'
        result[1].val == 'WORLD'
    }


    def 'should allow multiple use of the same process' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')

        MODULE.text = '''
            process producer {
                output: stdout
                shell: "echo Hello"
            }
            
            process consumer {
                input: file "foo"
                output: stdout
                shell:
                "cmd consumer 1"
            }
            
            process another_consumer {
                input: file "foo"
                output: stdout
                shell: "cmd consumer 2"
            }
            
            workflow flow1 {
                emit: producer | consumer | map { it.toUpperCase() }
            }
            
            workflow flow2 {
                emit: producer | another_consumer | map { it.toUpperCase() }
            }
            '''.stripIndent()

        when:
        def result = dsl_eval("""
            include { flow1; flow2 } from "$MODULE" 
  
            workflow { 
              flow1()
              flow2()
              emit: 
              flow1.out
              flow2.out
            }

            """)

        then:
        result[0].val == 'CMD CONSUMER 1'
        result[1].val == 'CMD CONSUMER 2'

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
            emit: foo.out
            emit: bar.out
        }
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result[0].val == 'HELLO'
        result[1].val == 'WORLD'
    }

    def 'should inherit module params' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        params.alpha = 'first'
        params.omega = 'last'
        
        process foo {
          output: val result
          exec:
            result = "$params.alpha $params.omega".toUpperCase()
        }     
        '''

        SCRIPT.text = """
        params.alpha = 'owner'
        include { foo } from "$MODULE" 

        workflow {
            foo()
            emit: foo.out
        }
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result.val == 'OWNER LAST'
    }

    def 'should override module params' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        params.alpha = 'first'
        params.omega = 'last'
        
        process foo {
          output: val result
          exec:
            result = "$params.alpha $params.omega".toUpperCase()
        }     
        '''

        SCRIPT.text = """
        params.alpha = 'owner'
        include { foo } from "$MODULE" params(alpha:'aaa', omega:'zzz')

        workflow {
            foo()
            emit: foo.out
        }
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result.val == 'AAA ZZZ'
    }

    def 'should extends module params' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        params.alpha = 'first'
        params.omega = 'last'
        
        process foo {
          output: val result
          exec:
            result = "$params.alpha $params.omega".toUpperCase()
        }     
        '''

        SCRIPT.text = """
        params.alpha = 'one' 
        params.omega = 'two'

        include { foo } from "$MODULE" addParams(omega:'zzz')

        workflow {
            foo()
            emit: foo.out
        }
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result.val == 'ONE ZZZ'
    }

    def 'should declare moduleDir path' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module/dir').createDirIfNotExists().resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = """
        def foo () { return true }
        assert moduleDir == file("$folder/module/dir")
        assert projectDir == file("$folder")
        assert launchDir == file('.')
        """

        SCRIPT.text = """
        include { foo } from "$MODULE" 

        assert moduleDir == file("$folder")
        assert projectDir == file("$folder")
        assert launchDir == file('.')
        
        workflow { true }
        """

        when:
        new MockScriptRunner()
                .setScript(SCRIPT)
                .execute()
        then:
        true
    }

    def 'should not allow unwrapped include' () {
        given:
        def folder = TestHelper.createInMemTempDir()
        def MODULE = folder.resolve('module.nf')
        def SCRIPT = folder.resolve('main.nf')

        MODULE.text = '''
        params.foo = 'x' 
        params.bar = 'y'
        
        process foo {
          output: stdout 
          script:
          /echo $params.foo $params.bar/
        }
        '''

        // inject params in the module
        // and invoke the process 'foo'
        SCRIPT.text = """
        include foo from "./module.nf" params(foo:'Hello', bar: 'world')
            
        workflow { 
            main: foo()
            emit: foo.out 
        }
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        def e = thrown(DeprecationException)
        e.message == "Unwrapped module inclusion is deprecated -- Replace `include foo from './MODULE/PATH'` with `include { foo } from './MODULE/PATH'`"

    }

    def 'should not allow include nested within a workflow' () {
        given:
        def folder = TestHelper.createInMemTempDir()
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
        new MockScriptRunner().setScript(SCRIPT).execute()
        then:
        def e = thrown(IllegalStateException)
        e.message == "Include statement is not allowed within a workflow definition"

    }

    def 'should should allow invoking function passing gstring' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')

        MODULE.text = '''
        def alpha(String str) {
            return str.reverse()
        }   
        '''

        when:
        def SCRIPT = """
        include { alpha } from "$MODULE"
        
        def x = "world"
        def y = "Hello \$x"
        
        return alpha(y)
        """
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        result == 'dlrow olleH'
    }
}
