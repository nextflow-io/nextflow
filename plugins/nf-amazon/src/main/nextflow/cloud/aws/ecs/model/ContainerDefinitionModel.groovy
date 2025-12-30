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
import software.amazon.awssdk.services.ecs.model.ContainerDefinition
import software.amazon.awssdk.services.ecs.model.KeyValuePair
import software.amazon.awssdk.services.ecs.model.LogConfiguration
import software.amazon.awssdk.services.ecs.model.LogDriver
import software.amazon.awssdk.services.ecs.model.ResourceRequirement
import software.amazon.awssdk.services.ecs.model.ResourceType

/**
 * Mutable wrapper model for building ECS ContainerDefinition.
 *
 * This class provides a builder-style interface for constructing container
 * definitions with Nextflow task requirements.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
@CompileStatic
class ContainerDefinitionModel {

    /**
     * Container name (always "main" for Nextflow tasks)
     */
    String name = 'main'

    /**
     * Container image (e.g., "ubuntu:latest")
     */
    String image

    /**
     * Command to execute (typically the bash wrapper script)
     */
    List<String> command

    /**
     * Entry point override
     */
    List<String> entryPoint

    /**
     * Working directory inside the container
     */
    String workingDirectory

    /**
     * Environment variables
     */
    Map<String, String> environment = [:]

    /**
     * Whether the container is essential (must be true for single-container tasks)
     */
    boolean essential = true

    /**
     * CloudWatch Logs group name
     */
    String logsGroup

    /**
     * CloudWatch Logs stream prefix
     */
    String logsStreamPrefix

    /**
     * AWS region for CloudWatch Logs
     */
    String logsRegion

    /**
     * CPU units for this container (soft limit)
     */
    Integer cpu

    /**
     * Memory in MiB for this container (hard limit)
     */
    Integer memory

    /**
     * Memory reservation in MiB (soft limit)
     */
    Integer memoryReservation

    /**
     * Build the ECS ContainerDefinition from this model.
     *
     * @param gpuCount Number of GPUs to allocate (null or 0 for none)
     * @return ContainerDefinition ready for task definition
     */
    ContainerDefinition toContainerDefinition(Integer gpuCount = null) {
        def builder = ContainerDefinition.builder()
            .name(name)
            .essential(essential)

        if (image)
            builder.image(image)

        if (command)
            builder.command(command)

        if (entryPoint)
            builder.entryPoint(entryPoint)

        if (workingDirectory)
            builder.workingDirectory(workingDirectory)

        // Add environment variables
        if (environment) {
            def envPairs = environment.collect { k, v ->
                KeyValuePair.builder().name(k).value(v).build()
            }
            builder.environment(envPairs)
        }

        // Configure CloudWatch Logs
        if (logsGroup) {
            def logOptions = [
                'awslogs-group': logsGroup,
                'awslogs-stream-prefix': logsStreamPrefix ?: 'nf'
            ]
            if (logsRegion)
                logOptions['awslogs-region'] = logsRegion

            builder.logConfiguration(
                LogConfiguration.builder()
                    .logDriver(LogDriver.AWSLOGS)
                    .options(logOptions)
                    .build()
            )
        }

        // Add CPU/memory limits if specified
        if (cpu)
            builder.cpu(cpu)

        if (memory)
            builder.memory(memory)

        if (memoryReservation)
            builder.memoryReservation(memoryReservation)

        // Add GPU resource requirements if specified
        if (gpuCount && gpuCount > 0) {
            builder.resourceRequirements(
                ResourceRequirement.builder()
                    .type(ResourceType.GPU)
                    .value(gpuCount.toString())
                    .build()
            )
        }

        return builder.build()
    }

    /**
     * Add an environment variable
     */
    ContainerDefinitionModel withEnvironment(String key, String value) {
        if (value != null)
            environment[key] = value
        return this
    }

    /**
     * Add multiple environment variables
     */
    ContainerDefinitionModel withEnvironment(Map<String, String> env) {
        if (env)
            environment.putAll(env)
        return this
    }

    /**
     * Configure CloudWatch Logs
     */
    ContainerDefinitionModel withLogging(String group, String prefix, String region) {
        this.logsGroup = group
        this.logsStreamPrefix = prefix
        this.logsRegion = region
        return this
    }
}
