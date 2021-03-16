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
                foo.emissionNext {
                    assert out[0] == 'echo hello'
                }
            }
            
            testflow test2 {
              when:
                bar(foo())
              then:
                bar.emissionNext {
                    assert out[0] == "say 'echo hello'"
                    assert meta.name == 'test2:bar'
                }
            }
        '''

        when:
        testflow_eval(SCRIPT)
        then:
        noExceptionThrown()

    }

    def 'should run test with tag' () {

        given:
        def SCRIPT = '''
            process foo {
                tag "$id"
                input: val(id)
                output: stdout
                script:
                "echo hello $id"
            }
        
            testflow test1 {
              when:
                foo( channel.of('A','B') )
              then:
                assert foo.outputsCount() == 1
                assert foo.emissionsCount() == 2
                
                foo.emissionWith(tag:'A') {
                    assert out[0] == 'echo hello A'
                }

                foo.emissionWith(tag:'B') {
                    assert out[0] == 'echo hello B'
                }
            }
        '''

        when:
        testflow_eval(SCRIPT)
        then:
        noExceptionThrown()

    }

    def 'should run test with index' () {

        given:
        def SCRIPT = '''
            process foo {
                tag "$id"
                input: val(id)
                output: stdout
                script:
                "echo hello $id"
            }
        
            testflow test1 {
              when:
                foo( channel.of('A','B','C') )
              then:
                assert foo.outputsCount() == 1
                assert foo.emissionsCount() == 3
                
                foo.emissionWith(index:1) {
                    assert out[0] == 'echo hello A'
                }

                foo.emissionWith(index:2) {
                    assert out[0] == 'echo hello B'
                }

                foo.emissionWith(index:3) {
                    assert out[0] == 'echo hello C'
                }
            }
        '''

        when:
        testflow_eval(SCRIPT)
        then:
        noExceptionThrown()

    }

}
