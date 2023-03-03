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

package io.seqera.wave.plugin.resolver

import io.seqera.wave.plugin.WaveClient
import nextflow.container.ContainerConfig
import nextflow.container.resolver.ContainerInfo
import nextflow.container.resolver.DefaultContainerResolver
import nextflow.executor.Executor
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WaveContainerResolverTest extends Specification {

    def 'should resolve docker container with wave container' () {
        given:
        def CONTAINER_NAME = "ubuntu:latest"
        def WAVE_CONTAINER = new ContainerInfo(CONTAINER_NAME, "wave.io/ubuntu:latest", "12345")
        def SINGULARITY_CONTAINER = new ContainerInfo('ubuntu:latest', '/some/singularity/ubuntu.img')
        and:
        def defaultResolver = Spy(DefaultContainerResolver)
        def executor = Mock(Executor)
        def resolver = Spy(new WaveContainerResolver(defaultResolver: defaultResolver))
        def task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) {
                getExecutor() >> executor
            }
        }

        when:
        def result = resolver.resolveImage(task, CONTAINER_NAME)
        then:
        resolver.client() >> Mock(WaveClient) { enabled()>>true }
        _ * task.getContainerConfig() >> Mock(ContainerConfig) { getEngine()>>'docker' }
        and:
        1 * resolver.waveContainer(task, CONTAINER_NAME) >> WAVE_CONTAINER
        and:
        result == WAVE_CONTAINER


        when:
        result = resolver.resolveImage(task, CONTAINER_NAME)
        then:
        resolver.client() >> Mock(WaveClient) { enabled()>>true }
        _ * task.getContainerConfig() >> Mock(ContainerConfig) { getEngine()>>'singularity' }
        and:
        1 * resolver.waveContainer(task, CONTAINER_NAME) >> WAVE_CONTAINER
        1 * defaultResolver.resolveImage(task, WAVE_CONTAINER.target) >> SINGULARITY_CONTAINER
        and:
        result == SINGULARITY_CONTAINER
    }

    def 'should validate container name' () {
        when:
        WaveContainerResolver.validateContainerRepo('ubuntu')
        then:
        noExceptionThrown()

        when:
        WaveContainerResolver.validateContainerRepo('ubuntu:latest')
        then:
        noExceptionThrown()

        when:
        WaveContainerResolver.validateContainerRepo('quay.io/wtsicgp/nanoseq:3.3.0')
        then:
        noExceptionThrown()

        when:
        WaveContainerResolver.validateContainerRepo('docker://quay.io/wtsicgp/nanoseq:3.3.0')
        then:
        thrown(IllegalArgumentException)
    }
}
