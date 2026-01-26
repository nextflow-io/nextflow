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

import static test.ScriptHelper.*

import java.nio.file.Files
import java.nio.file.Paths

import nextflow.NextflowMeta
import nextflow.SysEnv
import nextflow.extension.FilesEx
import nextflow.secret.SecretsLoader
import test.Dsl2Spec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BaseScriptTest extends Dsl2Spec {


    def 'should define implicit variables' () {

        given:
        def folder = Files.createTempDirectory('test')
        def WORK_DIR = folder.resolve('work')

        when:
        def config = [workDir: WORK_DIR]
        def script = folder.resolve('main.nf')
        script.text = '''
                result = [:]
                result.baseDir = baseDir
                result.projectDir = projectDir
                result.workDir = workDir 
                result.nextflow = nextflow
                result.workflow = workflow
                result.launchDir = launchDir 
                result.moduleDir = moduleDir
                result
                '''

        def result = runScript(script, config: config)

        then:
        result.baseDir == script.parent
        result.projectDir == script.parent
        result.workDir == WORK_DIR
        result.launchDir == Paths.get('.').toRealPath()
        // moduleDir uses real path (symlinks resolved) for consistent lookups
        result.moduleDir == script.parent
        result.workflow instanceof WorkflowMetadata
        result.nextflow == NextflowMeta.instance

        cleanup:
        script?.delete()
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

        when:
        script.text = '''
                secrets.FOO
                '''

        def result = runScript(script)

        then:
        result == 'ciao'

        cleanup:
        folder?.deleteDir()
        SysEnv.pop()
    }

}
