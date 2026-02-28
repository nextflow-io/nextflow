/*
 * Copyright 2013-2026, Seqera Labs
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

import test.Dsl2Spec

import static test.ScriptHelper.*

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
        def result = runScript(SCRIPT)

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
        def result = runScript(SCRIPT)

        then:
        result != null
        println "Path result: $result"
        println "Path result class: ${result?.getClass()}"

        cleanup:
        Files.deleteIfExists(tempFile)
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
        def result = runScript(SCRIPT)

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
        runScript(SCRIPT)

        then:
        def e = thrown(Exception)
        e.message.contains('Missing required parameter: --requiredParam')
    }
}
