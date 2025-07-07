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
 *
 */

package nextflow.k8s

import java.nio.file.Files

import nextflow.Session
import nextflow.executor.Executor
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class K8sWrapperBuilderTest extends Specification {

    def 'should render launcher script' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def sess = Mock(Session)
        def exec = Mock(Executor)
        def proc = Mock(TaskProcessor) { getSession() >> sess; getExecutor() >> exec }
        def config = new TaskConfig()
        def task = Mock(TaskRun) {
            getName() >> 'foo'
            getConfig() >> config
            getProcessor() >> proc
            getWorkDir() >> folder
            getInputFilesMap() >> [:]
            getOutputFilesNames() >> []
        }

        and:
        def builder = Spy(new K8sWrapperBuilder(task)) { getSecretsEnv() >> null; fixOwnership() >> false }

        when:
        def binding =  builder.makeBinding()

        then:
        binding.header_script == "NXF_CHDIR=${folder}"

        cleanup:
        folder?.deleteDir()
    }
}
