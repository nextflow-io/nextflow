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

package nextflow.executor.local

import nextflow.Session
import nextflow.SysEnv
import nextflow.container.DockerConfig
import nextflow.container.PodmanConfig
import spock.lang.Specification
import spock.lang.Unroll

class LocalExecutorTest extends Specification {

    @Unroll
    def 'should enable docker wrapper mode with config=#CONFIG windows=#WINDOWS env=#ENV fusion=#FUSION' () {
        given:
        SysEnv.push(ENV)
        def session = Mock(Session) { getContainerConfig() >> CONFIG }
        def executor = Spy(LocalExecutor) {
            isWindowsHost() >> WINDOWS
            isFusionEnabled() >> FUSION
        }
        executor.session = session

        expect:
        executor.isDockerWrapperEnabled() == EXPECTED
        executor.isContainerNative() == (EXPECTED || FUSION)

        cleanup:
        SysEnv.pop()

        where:
        CONFIG                              | WINDOWS | ENV                                         | FUSION | EXPECTED
        new DockerConfig([:])               | false   | [:]                                         | false  | false
        new DockerConfig(enabled: true)     | false   | [:]                                         | false  | false
        new DockerConfig(enabled: true)     | true    | [:]                                         | false  | true
        new DockerConfig(enabled: false)    | true    | [:]                                         | false  | false
        new DockerConfig(enabled: true)     | false   | [NXF_DOCKER_WRAPPER_IN_CONTAINER: 'true']    | false  | true
        new DockerConfig(enabled: true)     | true    | [:]                                         | true   | false
        new PodmanConfig(enabled: true)     | true    | [:]                                         | false  | false
        new PodmanConfig(enabled: true)     | false   | [NXF_DOCKER_WRAPPER_IN_CONTAINER: 'true']    | false  | false
    }

}
