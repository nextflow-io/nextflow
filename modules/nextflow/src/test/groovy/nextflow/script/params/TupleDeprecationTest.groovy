package nextflow.script.params

import spock.lang.Timeout
import test.Dsl2Spec
import test.MockScriptRunner

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class TupleDeprecationTest extends Dsl2Spec {

    def 'should not allow unqualified input file' () {
        given:
        def SCRIPT = '''
         
        process foo {
          input: 
          tuple 'x'
          /touch x/
        }
       
        workflow {
            foo()
        }
        '''

        when:
        new MockScriptRunner() .setScript(SCRIPT).execute()
        then:
        def e = thrown(DeprecationException)
        e.message == "Unqualified input file declaration has been deprecated - replace `tuple 'x',..` with `tuple path('x'),..`"
    }

    def 'should not allow unqualified input val' () {
        given:
        def SCRIPT = '''
         
        process foo {
          input: 
          tuple X
          /echo $X/
        }
       
        workflow {
            foo()
        }
        '''

        when:
        new MockScriptRunner() .setScript(SCRIPT).execute()
        then:
        def e = thrown(DeprecationException)
        e.message == "Unqualified input value declaration has been deprecated - replace `tuple X,..` with `tuple val(X),..`"
    }


    def 'should not allow unqualified output file' () {
        given:
        def SCRIPT = '''
         
        process foo {
          output: 
          tuple 'x'
          /touch x/
        }
       
        workflow {
            foo()
        }
        '''

        when:
        new MockScriptRunner() .setScript(SCRIPT).execute()
        then:
        def e = thrown(DeprecationException)
        e.message == "Unqualified output path declaration has been deprecated - replace `tuple 'x',..` with `tuple path('x'),..`"
    }

    def 'should not allow unqualified output value' () {
        given:
        def SCRIPT = '''
         
        process foo {
          output: 
          tuple X
          /echo hello/
        }
       
        workflow {
            foo()
        }
        '''

        when:
        new MockScriptRunner() .setScript(SCRIPT).execute()
        then:
        def e = thrown(DeprecationException)
        e.message == "Unqualified output value declaration has been deprecated - replace `tuple X,..` with `tuple val(X),..`"
    }
}
