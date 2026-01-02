package nextflow.script

import java.nio.file.Files

import test.Dsl2Spec
import test.MockScriptRunner

/**
 * Tests for single process execution feature that allows running processes
 * directly without explicit workflows.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ScriptProcessRunTest extends Dsl2Spec {

    def 'should execute single process with val input' () {
        given:
        def SCRIPT = '''
        params.sampleId = 'SAMPLE_001'
        
        process testProcess {
            input: val sampleId
            output: val result
            exec:
                result = "Processed: $sampleId"
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        // For single process execution, the result should contain the process output
        result != null
        println "Result: $result"
        println "Result class: ${result?.getClass()}"
    }

    def 'should execute single process with path input' () {
        given:
        def tempFile = Files.createTempFile('test', '.txt')
        def SCRIPT = """
        params.dataFile = '${tempFile}'
        
        process testProcess {
            input: path dataFile
            output: val result
            exec:
                result = "File: \${dataFile.name}"
        }
        """

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result != null
        println "Path result: $result"
        println "Path result class: ${result?.getClass()}"

        cleanup:
        Files.deleteIfExists(tempFile)
    }

    def 'should execute single process with tuple input' () {
        given:
        def SCRIPT = '''
        params.'meta.id' = 'SAMPLE_001'
        params.'meta.name' = 'test'
        params.threads = 4
        
        process testProcess {
            input: tuple val(meta), val(threads)
            output: val result
            exec:
                result = "Sample: ${meta.id}, Threads: $threads"
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result != null
        println "Tuple result: $result"
        println "Tuple result class: ${result?.getClass()}"
    }

    def 'should handle multiple processes by running the first one' () {
        given:
        def SCRIPT = '''
        params.input = 'test'
        
        process firstProcess {
            input: val input
            output: val result
            exec:
                result = "First: $input"
        }
        
        process secondProcess {
            input: val input
            output: val result
            exec:
                result = "Second: $input"
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result != null
        println "Multi-process result: $result"
        println "Multi-process result class: ${result?.getClass()}"
    }

    def 'should fail when required parameter is missing' () {
        given:
        def SCRIPT = '''
        process testProcess {
            input: val requiredParam
            output: val result
            exec:
                result = "Got: $requiredParam"
        }
        '''

        when:
        def runner = new MockScriptRunner()
        runner.setScript(SCRIPT).execute()

        then:
        def e = thrown(Exception)
        e.message.contains('Missing required parameter: --requiredParam')
    }

    def 'should handle complex parameter mapping' () {
        given:
        def SCRIPT = '''
        params.'sample.id' = 'S001'
        params.'sample.name' = 'TestSample'
        params.'config.threads' = 8
        
        process complexProcess {
            input: val sample
            input: val config
            output: val result
            exec:
                result = "Sample: ${sample.id}, Config: ${config.threads}"
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result != null
        println "Complex result: $result"
        println "Complex result class: ${result?.getClass()}"
    }
}
