/*
 * Copyright 2013-2024, Seqera Labs
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
import java.nio.file.Paths

import nextflow.NextflowMeta
import nextflow.Session
import nextflow.SysEnv
import nextflow.extension.FilesEx
import nextflow.secret.SecretsLoader
import test.Dsl2Spec
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BaseScriptTest extends Dsl2Spec {


    def 'should define implicit variables' () {

        given:
        def script = Files.createTempFile('test',null)
        and:
        def WORKFLOW = Mock(WorkflowMetadata)
        def WORK_DIR = Paths.get('/work/dir')
        def PROJECT_DIR = Paths.get('/some/base')
        and:
        def session = Mock(Session) {
            getBaseDir() >> PROJECT_DIR
            getWorkDir() >> WORK_DIR
            getWorkflowMetadata() >> WORKFLOW
        }
        def binding = new ScriptBinding([:])
        def parser = new ScriptParser(session)

        when:
        script.text = '''
                result = [:]
                result.baseDir = baseDir
                result.projectDir = projectDir
                result.workDir = workDir 
                result.nextflow = nextflow
                result.workflow = workflow
                result.launchDir = launchDir 
                result.moduleDir = moduleDir
                '''

        parser.setBinding(binding)
        parser.runScript(script)

        then:
        binding.result.baseDir == PROJECT_DIR
        binding.result.projectDir == PROJECT_DIR
        binding.result.workDir == WORK_DIR
        binding.result.launchDir == Paths.get('.').toRealPath()
        binding.result.moduleDir == script.parent
        binding.workflow == WORKFLOW
        binding.nextflow == NextflowMeta.instance

        cleanup:
        script?.delete()
    }

    def 'should use custom entry workflow' () {

        given:
        def script = Files.createTempFile('test',null)
        and:
        def session = Mock(Session)
        def binding = new ScriptBinding([:])
        def parser = new ScriptParser(session)

        when:
        script.text = '''
                workflow foo {
                }

                workflow {
                    error 'you were supposed to run foo!'
                }
                '''

        parser.setBinding(binding)
        parser.setEntryName('foo')
        parser.runScript(script)

        then:
        noExceptionThrown()

        cleanup:
        script?.delete()
    }

    def 'should use entry workflow from module' () {

        given:
        def folder = TestHelper.createInMemTempDir()
        def module = folder.resolve('module.nf')
        def script = folder.resolve('main.nf')
        and:
        def session = Mock(Session)
        def binding = new ScriptBinding([:])
        def parser = new ScriptParser(session)

        when:
        module.text = '''
                workflow foo {
                }
                '''

        script.text = '''
                include { foo } from './module.nf'

                workflow {
                    error 'you were supposed to run foo!'
                }
                '''

        parser.setBinding(binding)
        parser.setEntryName('foo')
        parser.runScript(script)

        then:
        noExceptionThrown()

        cleanup:
        folder?.delete()
    }

    def 'should resolve secret in a script' () {
        given:
        SecretsLoader.instance.reset()
        and:
        def folder = Files.createTempDirectory('test')
        def script = folder.resolve('main.nf')
        def secrets  = folder.resolve('store.json')
        and:
        secrets.text = '''
            [
              {
                "name": "FOO",
                "value": "ciao"
              }
            ]
            '''
        and:
        FilesEx.setPermissions(secrets, 'rw-------')
        SysEnv.push(NXF_SECRETS_FILE:secrets.toAbsolutePath().toString())
        and:
        def session = Mock(Session)
        def binding = new ScriptBinding([:])
        def parser = new ScriptParser(session)

        when:
        script.text = '''
                return secrets.FOO
                '''

        def result = parser
            .setBinding(binding)
            .runScript(script)
            .getResult()

        then:
        result == 'ciao'

        cleanup:
        folder?.deleteDir()
        SysEnv.pop()
    }

}
