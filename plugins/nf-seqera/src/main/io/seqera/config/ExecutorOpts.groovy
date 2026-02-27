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

package io.seqera.config

import groovy.transform.CompileStatic
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.script.dsl.Description
import nextflow.util.Duration

/**
 * Configuration for the Seqera executor.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Description("""
    The `seqera.executor` scope provides configuration for the Seqera compute executor.
""")
@CompileStatic
class ExecutorOpts implements ConfigScope {

    final RetryOpts retryPolicy

    @ConfigOption
    @Description("""
        The Seqera scheduler service endpoint URL.
    """)
    final String endpoint

    @ConfigOption
    @Description("""
        The AWS region for task execution (default: `eu-central-1`).
    """)
    final String region

    @ConfigOption
    @Description("""
        The EC2 key pair name for SSH access to instances.
    """)
    final String keyPairName

    @ConfigOption
    @Description("""
        The interval for batching task submissions (default: `1 sec`).
    """)
    final Duration batchFlushInterval

    @Description("""
        Machine/infrastructure requirements for session tasks.
    """)
    final MachineRequirementOpts machineRequirement

    @ConfigOption
    @Description("""
        Custom labels to apply to AWS resources for cost tracking and resource organization.
        Labels are propagated to ECS tasks, capacity providers, and EC2 instances.
    """)
    final Map<String, String> labels

    @ConfigOption
    @Description("""
        When `true`, automatically adds workflow metadata labels (e.g. project name,
        run name, session ID) with the `nextflow.io/` prefix to the session (default: `false`).
    """)
    final boolean autoLabels

    @ConfigOption
    @Description("""
        The resource prediction model to use for estimating task resource requirements
        based on historical execution metrics. Supported values: `qr/v1` (quantile regression).
        When not set, no resource estimation is applied.
    """)
    final String predictionModel

    @ConfigOption
    @Description("""
        Custom environment variables to apply to all tasks submitted by the Seqera executor.
        These are merged with the Fusion environment variables, with Fusion variables taking precedence.
    """)
    final Map<String, String> taskEnvironment

    /* required by config scope -- do not remove */

    ExecutorOpts() {}

    ExecutorOpts(Map opts) {
        this.retryPolicy = new RetryOpts(opts.retryPolicy as Map ?: Map.of())
        this.endpoint = opts.endpoint as String
        if (!endpoint)
            throw new IllegalArgumentException("Missing Seqera endpoint - make sure to specify 'seqera.executor.endpoint' settings")

        this.region = opts.region as String ?: "eu-central-1"
        this.keyPairName = opts.keyPairName as String
        this.batchFlushInterval = opts.batchFlushInterval
            ? Duration.of(opts.batchFlushInterval as String)
            : Duration.of('1 sec')
        // machine requirement settings
        this.machineRequirement = new MachineRequirementOpts(opts.machineRequirement as Map ?: Map.of())
        // labels for cost tracking
        this.labels = opts.labels as Map<String, String>
        this.autoLabels = opts.autoLabels as boolean ?: false
        // prediction model
        this.predictionModel = parsePredictionModel(opts.predictionModel as String)
        // custom task environment variables
        this.taskEnvironment = opts.taskEnvironment as Map<String, String>
    }

    private static final Set<String> VALID_PREDICTION_MODELS = Set.of('qr/v1')

    private static String parsePredictionModel(String value) {
        if( !value )
            return null
        if( !VALID_PREDICTION_MODELS.contains(value) )
            throw new IllegalArgumentException("Invalid prediction model '${value}'. Supported values: ${VALID_PREDICTION_MODELS.join(', ')}")
        return value
    }

    RetryOpts retryOpts() {
        this.retryPolicy
    }

    String getEndpoint() {
        return endpoint
    }

    String getRegion() {
        return region
    }

    String getKeyPairName() {
        return keyPairName
    }

    Duration getBatchFlushInterval() {
        return batchFlushInterval
    }

    MachineRequirementOpts getMachineRequirement() {
        return machineRequirement
    }

    Map<String, String> getLabels() {
        return labels
    }

    boolean getAutoLabels() {
        return autoLabels
    }

    String getPredictionModel() {
        return predictionModel
    }

    Map<String, String> getTaskEnvironment() {
        return taskEnvironment
    }
}
