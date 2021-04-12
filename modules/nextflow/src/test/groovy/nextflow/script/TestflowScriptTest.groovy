package nextflow.script

import groovy.xml.XmlSlurper
import test.Dsl2Spec
import test.MockScriptRunner

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TestflowScriptTest extends Dsl2Spec {

    protected ScriptRunner testflow_eval(String str, String fileName) {
        final runner = new MockScriptRunner()
        runner.setScript(str, fileName).testFlow()
        return runner
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
        def runner = testflow_eval(SCRIPT, "test_run.nf")

        then:
        noExceptionThrown()

        and: 'there is a XUnit report'
        def result = new XmlSlurper().parseText(runner.session.workDir.resolve("test-results/TEST-test_run.xml").text)
        result.@tests == 2
        result.@failures == 0
        result.@errors == 0
        result.testcase.size() == 2

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
        def runner = testflow_eval(SCRIPT, "test_run_with_tag.nf")

        then:
        noExceptionThrown()

        and: 'there is a XUnit report'
        def result = new XmlSlurper().parseText(runner.session.workDir.resolve("test-results/TEST-test_run_with_tag.xml").text)
        result.@tests == 1
        result.@failures == 0
        result.@errors == 0
        result.testcase.size() == 1

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
        def runner = testflow_eval(SCRIPT, "test_run_with_index.nf")

        then:
        noExceptionThrown()

        and: 'there is a XUnit report'
        def result = new XmlSlurper().parseText(runner.session.workDir.resolve("test-results/TEST-test_run_with_index.xml").text)
        result.@tests == 1
        result.@failures == 0
        result.@errors == 0
        result.testcase.size() == 1

    }

}
