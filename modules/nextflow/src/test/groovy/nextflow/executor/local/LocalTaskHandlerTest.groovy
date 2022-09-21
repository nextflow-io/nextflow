/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.executor.local

import java.nio.file.Path

import nextflow.Session
import nextflow.container.ContainerConfig
import nextflow.file.http.XPath
import nextflow.processor.TaskBean
import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LocalTaskHandlerTest extends Specification {

    def 'should create local process builder' () {
        given:
        def task = Mock(TaskRun) {
            getWorkDir() >> Path.of('/some/work/dir')
            getConfig() >> Mock(TaskConfig)
        }
        and:
        def handler = Spy(new LocalTaskHandler(task, Mock(LocalExecutor)))

        when:
        def builder = handler.createLaunchProcessBuilder()
        then:
        handler.fusionEnabled() >> false
        and:
        builder.command() == ['/bin/bash','-ue','.command.run']
        builder.directory() == new File('/some/work/dir')
        builder.redirectErrorStream()
        builder.redirectOutput().file() == new File('/some/work/dir/.command.log')
    }

    def 'should create fusion process builder' () {
        given:
        def WORK_DIR = XPath.get('http://some/work/dir')
        and:
        def session = Mock(Session) {
            getContainerConfig() >> new ContainerConfig([engine:'docker',enabled:true])
        }
        def bean = new TaskBean(workDir: WORK_DIR, inputFiles: [:])
        and:
        def task = Mock(TaskRun) {
            getContainer() >> 'ubuntu:latest'
            getWorkDir() >> WORK_DIR
            getConfig() >> Mock(TaskConfig)
            toTaskBean() >> bean
        }
        def executor = Mock(LocalExecutor) { getSession() >> session }
        and:
        def handler = Spy(new LocalTaskHandler(task, executor))

        when:
        def builder = handler.createLaunchProcessBuilder()
        then:
        handler.fusionEnabled() >> true
        and:
        builder.command() == ['sh','-c','docker run -i -e "NXF_FUSION_WORK=/fusion/http/some/work/dir" -e "NXF_FUSION_BUCKETS=http://some" --rm --privileged ubuntu:latest bash -o pipefail -c \'trap "{ ret=$?; cp .command.log /fusion/http/some/work/dir/.command.log||true; exit $ret; }" EXIT; bash /fusion/http/some/work/dir/.command.run 2>&1 | tee .command.log\'']
        builder.directory() == null
        builder.redirectErrorStream()
        builder.redirectOutput().file()

        cleanup:
        builder?.redirectOutput()?.file()?.delete()
    }
}
