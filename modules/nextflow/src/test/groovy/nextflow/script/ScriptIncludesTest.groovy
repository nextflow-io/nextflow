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

import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.NextflowMeta
import nextflow.exception.DuplicateModuleIncludeException
import spock.lang.Specification
import test.MockScriptRunner
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ScriptIncludesTest extends Specification {

    def setupSpec() { NextflowMeta.instance.enableDsl2() }
    def cleanupSpec() { NextflowMeta.instance.disableDsl2() }

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
        include "$MODULE" 
   
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
        
        workflow alpha(data) {
            foo(data)
            bar(foo.output)
        }
        
        '''

        SCRIPT.text = """
        include "$MODULE" 
   
        workflow {
            alpha('Hello')
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
        include "$MODULE"

        workflow alpha(data) {
            foo(data)
            bar(foo.output)
        }
   
        workflow {
            alpha('Hello')
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
        include "$MODULE"
   
        data = 'Hello'
        workflow {
            foo(data)
            bar(foo.output)
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
        include "$MODULE" 
   
        workflow {
            data = 'Hello'
            foo(data)
            bar(foo.output)
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
          output: stdout() 
          script:
          /echo Hello $sample/
        }        
        '''

        SCRIPT.text = """
        include "$MODULE" 
        hello_ch = Channel.from('world')
        
        workflow {
            foo(hello_ch)
        }
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        noExceptionThrown()
        result instanceof DataflowReadChannel
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
            set pairId, reads
          output: 
            stdout() 
          script:
            /echo sample=$sample pairId=$pairId reads=$reads/
        }
        '''

        SCRIPT.text = """
        include 'module.nf'

        workflow {
            ch1 = Channel.from('world')
            ch2 = Channel.value(['x', '/some/file'])
            foo(ch1, ch2)
        }
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        noExceptionThrown()
        result instanceof DataflowReadChannel
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
            gamma = 'world\'
            /nope/
        }
        
        process bar {
           input:
             val xx
             val yy 
           output:
             stdout()
           script:
            /echo $xx $yy/            
        }
        '''

        when:
        SCRIPT.text = """  
        include 'module.nf'        

        workflow {
            bar( foo('Ciao') )
        }
        """
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        noExceptionThrown()
        result instanceof DataflowReadChannel
        result.val == 'echo Ciao world'


        when:
        SCRIPT.text = """ 
        include 'module.nf'        
        
        workflow {
          (ch0, ch1) = foo('Ciao')
        }
        """
        runner = new MockScriptRunner()
        result = runner.setScript(SCRIPT).execute()
        then:
        noExceptionThrown()
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
          output: stdout() 
          script:
          /echo $params.foo $params.bar/
        }
        '''

        // inject params in the module
        // and invoke the process 'foo'
        SCRIPT.text = """     
        include "module.nf" params(foo:'Hello', bar: 'world')
            
        workflow { 
            foo() 
        }
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        noExceptionThrown()
        result instanceof DataflowReadChannel
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
        include 'module.nf'

        def str = foo('dlrow')
        return bar('Hello', str)
        """

        when:
        def runner = new TestScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
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
          output: stdout() 
          script:
          "echo $FOO"
        }
        '''

        SCRIPT.text = """ 
        include 'module.nf' params(x: 'Hola mundo')
        
        workflow {
            return foo()
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
        
        foo()
        '''

        SCRIPT.text = """ 
        include 'module.nf'
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
        include 'org/bio' 
        
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

    def 'should error on duplicate import' () {
        given:
        def folder = TestHelper.createInMemTempDir();
        def MOD1 = folder.resolve('mod1.nf')
        def MOD2 = folder.resolve('mod2.nf')
        def SCRIPT = folder.resolve('main.nf')

        MOD1.text = MOD2.text = '''     
        process foo {
            /hello/ 
        }      
        '''

        SCRIPT.text = """ 
        include 'mod1' 
        include 'mod2' 
        println 'x'
        """

        when:
        def runner = new MockScriptRunner()
        runner.setScript(SCRIPT).execute()
        then:
        def err = thrown(DuplicateModuleIncludeException)
        err.message == "A process with name 'foo' is already defined in the current context"
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
        include alpha from "$MODULE"
        alpha()
        """
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()
        then:
        result == 'this is alpha result'

        when:
        SCRIPT = """
        include alpha as FOO from "$MODULE"
        FOO()
        """
        runner = new MockScriptRunner()
        result = runner.setScript(SCRIPT).execute()
        then:
        result == 'this is alpha result'


        when:
        SCRIPT = """
        include alpha as FOO from "$MODULE"
        alpha()
        """
        runner = new MockScriptRunner()
        runner.setScript(SCRIPT).execute()
        then:
        thrown(MissingMethodException)


        when:
        SCRIPT = """
        include alpha from "$MODULE"
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
        include foo from "$MODULE" 
        include foo as bar from "$MODULE"  

        workflow {
            foo('Hello')
            bar('World')
            return [ foo.output, bar.output ]
        }
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result[0].val == 'HELLO'
        result[1].val == 'WORLD'
    }
}
