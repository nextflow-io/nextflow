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

import nextflow.Global
import nextflow.exception.AbortOperationException
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
        def result = runScript(SCRIPT, moduleRun: true)

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
        def result = runScript(SCRIPT, moduleRun: true)

        then:
        result != null
        println "Path result: $result"
        println "Path result class: ${result?.getClass()}"

        cleanup:
        Files.deleteIfExists(tempFile)
    }

    static final String ONE_PROCESS = '''
        process firstProcess {
            input: val input
            output: val result
            exec:
                result = "First: $input"
        }
        '''

    static final String TWO_PROCESSES = ONE_PROCESS + '''
        process secondProcess {
            input: val input
            output: val result
            exec:
                result = "Second: $input"
        }
        '''

    static final String PROCESS_AND_WORKFLOW = ONE_PROCESS + '''
        workflow testWorkflow {
            take: input
            main:
                firstProcess(input)
        }
        '''

    def 'should not execute a definition directly when it is not the only one, or without module run' () {
        when:
        runScript(SCRIPT, moduleRun: MODULE_RUN, config: [params: [input: 'test']])

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('No entry workflow specified')

        where:
        SCRIPT               | MODULE_RUN
        TWO_PROCESSES        | true
        PROCESS_AND_WORKFLOW | true
        // a single process is executable, but only via `nextflow module run`
        ONE_PROCESS          | false
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
        runScript(SCRIPT, moduleRun: true)

        then:
        def e = thrown(Exception)
        e.message.contains('Parameter `--requiredParam` is required but no value was provided')
    }

    def 'should cast boolean parameter to boolean' () {
        given:
        def folder = Files.createTempDirectory('test')

        def module = folder.resolve('module.nf')
        module.text = '''\
            process quant {
                input:
                val alignment_mode

                exec:
                assert alignment_mode instanceof Boolean
            }
            '''

        def spec = folder.resolve('meta.yml')
        spec.text = '''\
            input:
            - name: alignment_mode
              type: boolean
            '''

        when:
        runScript(module, moduleRun: true, config: [params: [alignment_mode: 'false']])
        then:
        Global.session.isSuccess()

        cleanup:
        folder?.deleteDir()
    }
}
