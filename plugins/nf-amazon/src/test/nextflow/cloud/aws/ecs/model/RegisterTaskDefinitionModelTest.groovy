/*
 * Copyright 2020-2024, Seqera Labs
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

package nextflow.cloud.aws.ecs.model

import software.amazon.awssdk.services.ecs.model.Compatibility
import software.amazon.awssdk.services.ecs.model.NetworkMode
import spock.lang.Specification

/**
 * Tests for {@link RegisterTaskDefinitionModel}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class RegisterTaskDefinitionModelTest extends Specification {

    def 'should create model with defaults'() {
        when:
        def model = RegisterTaskDefinitionModel.create()

        then:
        model.containerDefinitions.size() == 1
        model.containerDefinitions.first().name == 'main'
    }

    def 'should compute cache key'() {
        given:
        def model = RegisterTaskDefinitionModel.create()
        model.cpu = '1024'
        model.memory = '2048'
        model.containerDefinitions.first().image = 'ubuntu:latest'

        when:
        def key1 = model.computeCacheKey()

        then:
        key1 == 'ubuntu:latest:1024:2048:0:30'

        when:
        model.gpuCount = 2
        model.ephemeralStorageGiB = 100
        def key2 = model.computeCacheKey()

        then:
        key2 == 'ubuntu:latest:1024:2048:2:100'
    }

    def 'should build request'() {
        given:
        def model = RegisterTaskDefinitionModel.create()
        model.family = 'nf-test'
        model.cpu = '1024'
        model.memory = '2048'
        model.executionRoleArn = 'arn:aws:iam::123:role/exec'
        model.taskRoleArn = 'arn:aws:iam::123:role/task'
        model.containerDefinitions.first().image = 'ubuntu:latest'

        when:
        def request = model.toRequest()

        then:
        request.family() == 'nf-test'
        request.cpu() == '1024'
        request.memory() == '2048'
        request.executionRoleArn() == 'arn:aws:iam::123:role/exec'
        request.taskRoleArn() == 'arn:aws:iam::123:role/task'
        request.networkMode() == NetworkMode.AWSVPC
        request.requiresCompatibilities() == [Compatibility.EC2]
        request.containerDefinitions().size() == 1
        request.containerDefinitions().first().name() == 'main'
        request.containerDefinitions().first().image() == 'ubuntu:latest'
    }

    def 'should build request with ephemeral storage'() {
        given:
        def model = RegisterTaskDefinitionModel.create()
        model.family = 'nf-test'
        model.cpu = '1024'
        model.memory = '2048'
        model.ephemeralStorageGiB = 100
        model.containerDefinitions.first().image = 'ubuntu:latest'

        when:
        def request = model.toRequest()

        then:
        request.ephemeralStorage() != null
        request.ephemeralStorage().sizeInGiB() == 100
    }

    def 'should not include ephemeral storage if 30 or less'() {
        given:
        def model = RegisterTaskDefinitionModel.create()
        model.family = 'nf-test'
        model.cpu = '1024'
        model.memory = '2048'
        model.ephemeralStorageGiB = 30
        model.containerDefinitions.first().image = 'ubuntu:latest'

        when:
        def request = model.toRequest()

        then:
        request.ephemeralStorage() == null
    }
}
