package nextflow.file.http

import com.github.tomakehurst.wiremock.junit.WireMockRule
import com.github.tomjankes.wiremock.WireMockGroovy
import nextflow.NextflowMeta
import org.junit.Rule
import test.Dsl2Spec
import test.MockScriptRunner
import test.OutputCapture

import java.nio.file.Files
import java.nio.file.Path


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
class FunctionalTest extends Dsl2Spec{

    @Rule
    OutputCapture capture = new OutputCapture()

    @Rule
    WireMockRule wireMockRule = new WireMockRule(18080)

    Path workDir

    def setup(){
        workDir = Files.createTempDirectory("nf")
    }

    def cleanup(){
        workDir.deleteDir()
    }

    /**
     * test passing values through environment variables
     */
    def 'test query params' () {

        given:
        def wireMock = new WireMockGroovy(18080)
        wireMock.stub {
            request {
                method "GET"
                url "/index.html?format=json"
            }
            response {
                status 200
                body '{ "test": "ok"}'
            }
        }
        wireMock.stub {
            request {
                method "GET"
                url "/index.html"
            }
            response {
                status 200
                body '{ "test": "error"}'
            }
        }

        def config = [process:[executor:'local'], workDir:workDir]

        def text = '''
            process hello {
                debug            
                input:
                    path finput            
                output:
                    stdout
                script:                            
                    println finput.target.toFile().text                 
            }
            
            workflow {
                main:
                    def json1 = file("http://localhost:18080/index.html?format=json")
                    def ch = Channel.fromPath(json1, glob:false) 
                    hello(ch)
                emit:
                    hello.out 
            }
        '''
        when:
        def runner = new MockScriptRunner(config)
        def result = runner.setScript(text).execute()
        def stdout = capture.toString()
        then:
        stdout.indexOf('{ "test": "ok"}') != -1
    }

}
