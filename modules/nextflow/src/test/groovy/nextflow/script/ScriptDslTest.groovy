package nextflow.script

import spock.lang.Specification

import nextflow.NextflowMeta
import org.codehaus.groovy.control.MultipleCompilationErrorsException
import test.MockScriptRunner
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ScriptDslTest extends Specification {

    def setupSpec() { NextflowMeta.instance.enableDsl2() }
    def cleanupSpec() { NextflowMeta.instance.disableDsl2() }

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
        
        workflow alpha(data) {
            foo(data)
            bar(foo.output)
        }
   
        workflow {
           alpha('Hello')
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result.val == 'HELLO MUNDO'
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
           hello_out = hello()
           hello_out.map { it.toUpperCase()  }
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result.val == 'HELLO'
    }

}
