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


import groovy.transform.CompileStatic
import software.amazon.awssdk.services.batch.model.JobDefinitionType
import software.amazon.awssdk.services.batch.model.PlatformCapability
import software.amazon.awssdk.services.batch.model.RegisterJobDefinitionRequest

/**
 * Custom mutable RegisterJobDefinitionRequest class that allows subclasses to modify the request
 * before converting it to the immutable AWS SDK object.
 *
 * This is a mutable version of {@link RegisterJobDefinitionRequest} required
 * to simplify the extension of container settings in the AWS Batch executor
 * and its sub-classes (e.g. nf-xpack).
 */
@CompileStatic
class RegisterJobDefinitionModel {

    private String jobDefinitionName

    private JobDefinitionType type

    private List<PlatformCapability> platformCapabilities

    private ContainerPropertiesModel containerProperties

    private Map<String,String> parameters

    private Map<String,String> tags

    RegisterJobDefinitionModel jobDefinitionName(String value) {
        this.jobDefinitionName = value
        return this
    }

    RegisterJobDefinitionModel type(JobDefinitionType value) {
        this.type = value
        return this
    }

    RegisterJobDefinitionModel platformCapabilities(List<PlatformCapability> value) {
        this.platformCapabilities = value
        return this
    }

    RegisterJobDefinitionModel containerProperties(ContainerPropertiesModel value) {
        this.containerProperties = value
        return this
    }

    RegisterJobDefinitionModel parameters(Map<String,String> value) {
        this.parameters = value
        return this
    }

    RegisterJobDefinitionModel tags(Map<String,String> value) {
        this.tags = value
        return this
    }

    RegisterJobDefinitionModel addTagsEntry(String key, String value) {
        if( this.tags==null )
            this.tags = new LinkedHashMap<>()
        this.tags.put(key, value)
        return this
    }

    String getJobDefinitionName() {
        return jobDefinitionName
    }

    JobDefinitionType getType() {
        return type
    }

    List<PlatformCapability> getPlatformCapabilities() {
        return platformCapabilities
    }

    ContainerPropertiesModel getContainerProperties() {
        return containerProperties
    }

    Map<String, String> getParameters() {
        return parameters
    }

    Map<String, String> getTags() {
        return tags
    }

    RegisterJobDefinitionRequest toBatchRequest() {
        final builder = RegisterJobDefinitionRequest.builder()

        if (jobDefinitionName)
            builder.jobDefinitionName(jobDefinitionName)
        if (type)
            builder.type(type)
        if (platformCapabilities)
            builder.platformCapabilities(platformCapabilities)
        if (containerProperties)
            builder.containerProperties(containerProperties.toBatchContainerProperties())
        if (parameters)
            builder.parameters(parameters)
        if (tags)
            builder.tags(tags)

        return (RegisterJobDefinitionRequest) builder.build()
    }

    @Override
    String toString() {
        return "RegisterJobDefinitionModel{" +
            "jobDefinitionName='" + jobDefinitionName + '\'' +
            ", type=" + type +
            ", platformCapabilities=" + platformCapabilities +
            ", containerProperties=" + containerProperties +
            ", parameters=" + parameters +
            ", tags=" + tags +
            '}';
    }
}
