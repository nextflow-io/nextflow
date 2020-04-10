/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BaseScriptTest extends Specification {


    def 'should define implicit variables' () {

        given:
        def script = Files.createTempFile('test',null)

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
        binding.result.baseDir ==PROJECT_DIR
        binding.result.projectDir == PROJECT_DIR
        binding.result.workDir == WORK_DIR
        binding.result.launchDir == Paths.get('.').toRealPath()
        binding.result.moduleDir == script.parent
        binding.workflow == WORKFLOW
        binding.nextflow == NextflowMeta.instance

        cleanup:
        script?.delete()
    }


}
