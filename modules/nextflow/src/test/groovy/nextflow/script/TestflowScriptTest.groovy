package nextflow.script

import test.Dsl2Spec
import test.MockScriptRunner

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TestflowScriptTest extends Dsl2Spec {

    protected testflow_eval(String str) {
        new MockScriptRunner() .setScript(str).testFlow()
    }

    def 'should run test' () {

        given:
        def SCRIPT = '''
            process foo {
                output: stdout
                script:
                'echo hello'
            }
            
            process bar {
                input: val x 
                output: stdout 
                script:
                "say '$x'"
            }
        
            testflow test1 {
              when:
                foo()
              then:
                assert foo.outputsCount() == 1
                assert foo.emissionsCount() == 1
                foo.withEmission(0) {
                    assert out[0] == 'echo hello'
                }
            }
            
            testflow test2 {
              when:
                bar(foo())
              then:
                bar.withEmission(0) {
                    assert out[0] == "say 'echo hello'"
                }
            }
        '''

        when:
        testflow_eval(SCRIPT)
        then:
        noExceptionThrown()

    }

}
