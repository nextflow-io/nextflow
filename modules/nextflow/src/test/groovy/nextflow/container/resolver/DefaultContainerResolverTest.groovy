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

package nextflow.container.resolver

import nextflow.container.ContainerConfig
import nextflow.executor.Executor
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DefaultContainerResolverTest extends Specification {

    def 'should return empty for null container' () {
        given:
        def resolver = new DefaultContainerResolver()
        expect:
        resolver.resolveImage(Mock(TaskRun), null) == ContainerInfo.EMPTY
    }

    def 'should resolve image' () {
        given:
        def resolver = new DefaultContainerResolver()
        def task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) { getExecutor() >> Mock(Executor)}
        }
        and:

        when:
        def result = resolver.resolveImage(task, 'ubuntu:latest')
        then:
        1 * task.getContainerConfig() >> new ContainerConfig([engine:'docker', enabled:true, registry:'quay.io'])
        and:
        result.source == 'ubuntu:latest'
        result.target == 'quay.io/ubuntu:latest'
        result.hashKey == 'quay.io/ubuntu:latest'
    }


}
