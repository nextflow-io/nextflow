package nextflow.script


import nextflow.NextflowMeta
import org.codehaus.groovy.control.MultipleCompilationErrorsException
import spock.lang.Specification
import spock.lang.Timeout
import test.MockScriptRunner
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class ScriptDslTest extends Specification {

    def setupSpec() { NextflowMeta.instance.enableDsl2() }
    def cleanupSpec() { NextflowMeta.instance.disableDsl2() }

    def 'should execute basic workflow' () {
        given:
        def SCRIPT = '''
   
        workflow {
            emit: result 
            main:
            result = 'Hello world'
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result[0].val == 'Hello world'

    }

    def 'should execute emit' () {
        given:
        def SCRIPT = '''
   
        workflow {
            emit:
            result = 'Hello world'
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result[0].val == 'Hello world'

    }

    def 'should emit expression' () {
        given:
        def SCRIPT = '''
         
        def foo() { 'Hello world' } 
       
        workflow {
            emit:
            foo().toUpperCase()
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result[0].val == 'HELLO WORLD'

    }

    def 'should emit process out' () {
        given:
        def SCRIPT = '''
         
        process foo {
          output: val x 
          exec: x = 'Hello'
        }
       
        workflow {
            main: foo()
            emit: foo.out
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result[0].val == 'Hello'

    }



    def 'should define processes and workflow' () {
        given:
        def SCRIPT = '''
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
            get: 
                data
            
            main:
                foo(data)
                bar(foo.out)
                
            emit: 
                x = bar.out 
            
        }
   
        workflow {
            main: alpha('Hello') 
            emit: x = alpha.out 
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result[0].val == 'HELLO MUNDO'
    }



    def 'should access nextflow enabling property' () {
        given:
        def SCRIPT = '''
        return nextflow.preview.dsl 
        '''
        when:
        def result = new MockScriptRunner().setScript(SCRIPT).execute()
        then:
        result == 2
    }


    def 'should not allow function with reserved identifier' () {

        given:
        def SCRIPT = """ 
            def main() { println 'ciao' }
        """

        when:
        new MockScriptRunner()
                .setScript(SCRIPT)
                .execute()
        then:
        def err = thrown(MultipleCompilationErrorsException)
        err.message.contains('Identifier `main` is reserved for internal use')
    }

    def 'should not allow process with reserved identifier' () {

        given:
        def SCRIPT = """ 
            process main {
              /echo ciao/
            }
        """

        when:
        new MockScriptRunner()
                .setScript(SCRIPT)
                .execute()
        then:
        def err = thrown(MultipleCompilationErrorsException)
        err.message.contains('Identifier `main` is reserved for internal use')
    }

    def 'should not allow workflow with reserved identifier' () {

        given:
        def SCRIPT = """ 
            workflow main {
              /echo ciao/
            }
        """

        when:
        new MockScriptRunner()
                .setScript(SCRIPT)
                .execute()
        then:
        def err = thrown(MultipleCompilationErrorsException)
        err.message.contains('Identifier `main` is reserved for internal use')
    }

    def 'should not allow duplicate workflow keyword' () {
        given:
        def SCRIPT = """ 
            workflow {
              /echo ciao/
            }
            
            workflow {
              /echo miao/
            }
        """

        when:
        new MockScriptRunner()
                .setScript(SCRIPT)
                .execute()
        then:
        def err = thrown(MultipleCompilationErrorsException)
        err.message.contains('Duplicate entry workflow definition')
    }

    def 'should apply operator to process result' () {
        given:
        def SCRIPT = '''
        process hello {
          output: val result
          exec:
            result = "Hello"
        }     
        
        workflow {
           main: hello()
           emit: hello.out.map { it.toUpperCase()  }
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result[0].val == 'HELLO'
    }

}
