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

package nextflow.cloud.aws.batch.model

import software.amazon.awssdk.services.batch.model.JobDefinitionType
import software.amazon.awssdk.services.batch.model.PlatformCapability
import software.amazon.awssdk.services.batch.model.RegisterJobDefinitionRequest
import spock.lang.Specification

/**
 * @author Nextflow Authors
 */
class RegisterJobDefinitionModelTest extends Specification {

    def 'should create empty model'() {
        when:
        def model = new RegisterJobDefinitionModel()

        then:
        model.jobDefinitionName == null
        model.type == null
        model.platformCapabilities == null
        model.containerProperties == null
        model.parameters == null
        model.tags == null
    }

    def 'should set and get job definition name'() {
        given:
        def model = new RegisterJobDefinitionModel()

        when:
        def result = model.jobDefinitionName('test-job-def')

        then:
        result == model
        model.jobDefinitionName == 'test-job-def'
    }

    def 'should set and get type'() {
        given:
        def model = new RegisterJobDefinitionModel()

        when:
        def result = model.type(JobDefinitionType.CONTAINER)

        then:
        result == model
        model.type == JobDefinitionType.CONTAINER
    }

    def 'should set and get platform capabilities'() {
        given:
        def model = new RegisterJobDefinitionModel()
        def capabilities = [PlatformCapability.EC2, PlatformCapability.FARGATE]

        when:
        def result = model.platformCapabilities(capabilities)

        then:
        result == model
        model.platformCapabilities == capabilities
        model.platformCapabilities.size() == 2
        model.platformCapabilities.contains(PlatformCapability.EC2)
        model.platformCapabilities.contains(PlatformCapability.FARGATE)
    }

    def 'should set and get container properties'() {
        given:
        def model = new RegisterJobDefinitionModel()
        def containerProps = new ContainerPropertiesModel()

        when:
        def result = model.containerProperties(containerProps)

        then:
        result == model
        model.containerProperties == containerProps
    }

    def 'should set and get parameters'() {
        given:
        def model = new RegisterJobDefinitionModel()
        def params = ['key1': 'value1', 'key2': 'value2']

        when:
        def result = model.parameters(params)

        then:
        result == model
        model.parameters == params
        model.parameters.size() == 2
        model.parameters['key1'] == 'value1'
        model.parameters['key2'] == 'value2'
    }

    def 'should set and get tags'() {
        given:
        def model = new RegisterJobDefinitionModel()
        def tags = ['env': 'test', 'project': 'nextflow']

        when:
        def result = model.tags(tags)

        then:
        result == model
        model.tags == tags
        model.tags.size() == 2
        model.tags['env'] == 'test'
        model.tags['project'] == 'nextflow'
    }

    def 'should add tag entry when tags is null'() {
        given:
        def model = new RegisterJobDefinitionModel()

        when:
        def result = model.addTagsEntry('key1', 'value1')

        then:
        result == model
        model.tags != null
        model.tags.size() == 1
        model.tags['key1'] == 'value1'
        model.tags instanceof LinkedHashMap
    }

    def 'should add tag entry when tags already exists'() {
        given:
        def model = new RegisterJobDefinitionModel()
        model.tags(['existing': 'tag'])

        when:
        def result = model.addTagsEntry('new', 'value')

        then:
        result == model
        model.tags.size() == 2
        model.tags['existing'] == 'tag'
        model.tags['new'] == 'value'
    }

    def 'should handle multiple tag entries'() {
        given:
        def model = new RegisterJobDefinitionModel()

        when:
        model.addTagsEntry('key1', 'value1')
             .addTagsEntry('key2', 'value2')
             .addTagsEntry('key3', 'value3')

        then:
        model.tags.size() == 3
        model.tags['key1'] == 'value1'
        model.tags['key2'] == 'value2'
        model.tags['key3'] == 'value3'
    }

    def 'should handle tag entry overwrite'() {
        given:
        def model = new RegisterJobDefinitionModel()

        when:
        model.addTagsEntry('key1', 'value1')
             .addTagsEntry('key1', 'value2')

        then:
        model.tags.size() == 1
        model.tags['key1'] == 'value2'
    }

    def 'should support method chaining'() {
        given:
        def model = new RegisterJobDefinitionModel()
        def containerProps = new ContainerPropertiesModel()
        def capabilities = [PlatformCapability.EC2]
        def params = ['param1': 'value1']
        def tags = ['tag1': 'value1']

        when:
        def result = model
            .jobDefinitionName('test-job')
            .type(JobDefinitionType.CONTAINER)
            .platformCapabilities(capabilities)
            .containerProperties(containerProps)
            .parameters(params)
            .tags(tags)
            .addTagsEntry('tag2', 'value2')

        then:
        result == model
        model.jobDefinitionName == 'test-job'
        model.type == JobDefinitionType.CONTAINER
        model.platformCapabilities == capabilities
        model.containerProperties == containerProps
        model.parameters == params
        model.tags.size() == 2
        model.tags['tag1'] == 'value1'
        model.tags['tag2'] == 'value2'
    }

    def 'should handle empty collections'() {
        given:
        def model = new RegisterJobDefinitionModel()

        when:
        model.platformCapabilities([])
             .parameters([:])
             .tags([:])

        then:
        model.platformCapabilities == []
        model.parameters == [:]
        model.tags == [:]
    }

    def 'should convert to RegisterJobDefinitionRequest with all fields'() {
        given:
        def model = new RegisterJobDefinitionModel()
        def containerProps = new ContainerPropertiesModel()
        containerProps.image('ubuntu:20.04')
        def capabilities = [PlatformCapability.EC2, PlatformCapability.FARGATE]
        def params = ['param1': 'value1', 'param2': 'value2']
        def tags = ['tag1': 'value1', 'tag2': 'value2']

        when:
        model.jobDefinitionName('test-job-def')
             .type(JobDefinitionType.CONTAINER)
             .platformCapabilities(capabilities)
             .containerProperties(containerProps)
             .parameters(params)
             .tags(tags)

        def request = model.toBatchRequest()

        then:
        request instanceof RegisterJobDefinitionRequest
        request.jobDefinitionName() == 'test-job-def'
        request.type() == JobDefinitionType.CONTAINER
        request.platformCapabilities() == capabilities
        request.containerProperties() != null
        request.containerProperties().image() == 'ubuntu:20.04'
        request.parameters() == params
        request.tags() == tags
    }

    def 'should convert to RegisterJobDefinitionRequest with null fields'() {
        given:
        def model = new RegisterJobDefinitionModel()

        when:
        def request = model.toBatchRequest()

        then:
        request instanceof RegisterJobDefinitionRequest
        !request.jobDefinitionName()
        !request.type()
        !request.platformCapabilities()
        !request.containerProperties()
        !request.parameters()
        !request.tags()
    }

    def 'should convert to RegisterJobDefinitionRequest with minimal fields'() {
        given:
        def model = new RegisterJobDefinitionModel()
        def containerProps = new ContainerPropertiesModel()
        containerProps.image('nginx')

        when:
        model.jobDefinitionName('minimal-job')
             .type(JobDefinitionType.CONTAINER)
             .containerProperties(containerProps)

        def request = model.toBatchRequest()

        then:
        request instanceof RegisterJobDefinitionRequest
        request.jobDefinitionName() == 'minimal-job'
        request.type() == JobDefinitionType.CONTAINER
        request.containerProperties() != null
        request.containerProperties().image() == 'nginx'
        !request.platformCapabilities()
        !request.parameters()
        !request.tags()
    }

    def 'should convert to RegisterJobDefinitionRequest with empty collections'() {
        given:
        def model = new RegisterJobDefinitionModel()
        def containerProps = new ContainerPropertiesModel()

        when:
        model.jobDefinitionName('empty-collections-job')
             .type(JobDefinitionType.CONTAINER)
             .platformCapabilities([])
             .containerProperties(containerProps)
             .parameters([:])
             .tags([:])

        def request = model.toBatchRequest()

        then:
        request instanceof RegisterJobDefinitionRequest
        request.jobDefinitionName() == 'empty-collections-job'
        request.type() == JobDefinitionType.CONTAINER
        request.platformCapabilities() == []
        request.containerProperties() != null
        request.parameters() == [:]
        request.tags() == [:]
    }

    def 'should handle chaining with toBatchRequest'() {
        given:
        def model = new RegisterJobDefinitionModel()
        def containerProps = new ContainerPropertiesModel()
        containerProps.image('alpine')

        when:
        def request = model
            .jobDefinitionName('chained-job')
            .type(JobDefinitionType.CONTAINER)
            .containerProperties(containerProps)
            .addTagsEntry('env', 'test')
            .addTagsEntry('project', 'nextflow')
            .toBatchRequest()

        then:
        request instanceof RegisterJobDefinitionRequest
        request.jobDefinitionName() == 'chained-job'
        request.type() == JobDefinitionType.CONTAINER
        request.containerProperties().image() == 'alpine'
        request.tags().size() == 2
        request.tags()['env'] == 'test'
        request.tags()['project'] == 'nextflow'
    }
}
