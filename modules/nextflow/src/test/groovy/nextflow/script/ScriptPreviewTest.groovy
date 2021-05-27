package nextflow.script


import nextflow.NextflowMeta
import spock.lang.Specification
import test.MockScriptRunner
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ScriptPreviewTest extends Specification {

    def setupSpec() { NextflowMeta.instance.enableDsl2(true) }

    def cleanupSpec() { NextflowMeta.instance.disableDsl2() }


    def 'should allow unwrapped include' () {
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
        noExceptionThrown()
        result.val == 'echo Hello world'

    }

}
