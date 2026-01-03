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

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import software.amazon.awssdk.services.ecs.model.Compatibility
import software.amazon.awssdk.services.ecs.model.EphemeralStorage
import software.amazon.awssdk.services.ecs.model.NetworkMode
import software.amazon.awssdk.services.ecs.model.RegisterTaskDefinitionRequest

/**
 * Mutable wrapper model for building ECS RegisterTaskDefinition requests.
 *
 * This class provides a builder-style interface for constructing task definitions
 * with Nextflow process requirements mapped to ECS task definition fields.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
@CompileStatic
class RegisterTaskDefinitionModel {

    /**
     * Task definition family name (e.g., "nf-ubuntu-latest")
     */
    String family

    /**
     * Task execution role ARN - used by ECS agent to pull images and write logs
     */
    String executionRoleArn

    /**
     * Task role ARN - used by the container for AWS API calls (e.g., S3 access)
     */
    String taskRoleArn

    /**
     * CPU units (1024 = 1 vCPU)
     */
    String cpu

    /**
     * Memory in MiB
     */
    String memory

    /**
     * Container definitions (typically just one "main" container)
     */
    List<ContainerDefinitionModel> containerDefinitions = []

    /**
     * Ephemeral storage size in GiB (30-200 for Fargate, 30-16384 for Managed Instances)
     */
    Integer ephemeralStorageGiB

    /**
     * Number of GPUs required (for EC2/Managed Instances only)
     */
    Integer gpuCount

    /**
     * Build the ECS RegisterTaskDefinitionRequest from this model.
     *
     * @return RegisterTaskDefinitionRequest ready for AWS SDK call
     */
    RegisterTaskDefinitionRequest toRequest() {
        def builder = RegisterTaskDefinitionRequest.builder()
            .family(family)
            .networkMode(NetworkMode.AWSVPC)
            .requiresCompatibilities(Compatibility.EC2)  // Managed Instances uses EC2 compatibility

        if (executionRoleArn)
            builder.executionRoleArn(executionRoleArn)

        if (taskRoleArn)
            builder.taskRoleArn(taskRoleArn)

        if (cpu)
            builder.cpu(cpu)

        if (memory)
            builder.memory(memory)

        // Build container definitions
        if (containerDefinitions) {
            def containers = containerDefinitions.collect { it.toContainerDefinition(gpuCount) }
            builder.containerDefinitions(containers)
        }

        // Configure ephemeral storage for Managed Instances
        if (ephemeralStorageGiB && ephemeralStorageGiB > 30) {
            builder.ephemeralStorage(
                EphemeralStorage.builder()
                    .sizeInGiB(ephemeralStorageGiB)
                    .build()
            )
        }

        return builder.build() as RegisterTaskDefinitionRequest
    }

    /**
     * Compute a cache key for this task definition based on resource requirements.
     * Used to avoid re-registering identical task definitions.
     *
     * @return Cache key string
     */
    String computeCacheKey() {
        def container = containerDefinitions ? containerDefinitions.first() : null
        def image = container?.image ?: 'unknown'
        def gpu = gpuCount ?: 0
        def disk = ephemeralStorageGiB ?: 30
        return "${image}:${cpu}:${memory}:${gpu}:${disk}"
    }

    /**
     * Create a model with default container named "main"
     */
    static RegisterTaskDefinitionModel create() {
        def model = new RegisterTaskDefinitionModel()
        model.containerDefinitions = [new ContainerDefinitionModel(name: 'main')]
        return model
    }
}
