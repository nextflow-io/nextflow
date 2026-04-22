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

    static final Set<String> VALID_AUTO_LABELS = Collections.unmodifiableSet(new LinkedHashSet<>([
        'projectName', 'userName', 'runName', 'sessionId', 'resume',
        'revision', 'commitId', 'repository', 'manifestName',
        'runtimeVersion', 'workflowId', 'workspaceId', 'computeEnvId'
    ]))

    final RetryOpts retryPolicy

    @ConfigOption
    @Description("""
        The Seqera scheduler service endpoint URL.
    """)
    final String endpoint

    @ConfigOption
    @Description("""
        The compute backend provider type (e.g. `aws`, `local`).
        When specified, used together with region to select the matching compute environment.
    """)
    final String provider

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
        Automatically attach workflow metadata labels (with the `nextflow.io/` and
        `seqera.io/platform/` prefixes) to the session. Accepts:
          - `true`: include all available metadata labels
          - `false` (default): disable
          - a list or comma-separated string of short names: e.g.
            `['runName', 'projectName']` or `'runName,projectName'`
        Valid names: `projectName`, `userName`, `runName`, `sessionId`, `resume`,
        `revision`, `commitId`, `repository`, `manifestName`, `runtimeVersion`,
        `workflowId`, `workspaceId`, `computeEnvId`.
    """)
    final Set<String> autoLabels

    @ConfigOption
    @Description("""
        The resource prediction model to use for estimating task resource requirements
        based on historical execution metrics. Supported values: `qr/v1`, `qr/v2` (quantile regression).
        When not set, no resource estimation is applied.
    """)
    final String predictionModel

    @ConfigOption
    @Description("""
        Custom environment variables to apply to all tasks submitted by the Seqera executor.
        These are merged with the Fusion environment variables, with Fusion variables taking precedence.
    """)
    final Map<String, String> taskEnvironment

    @ConfigOption
    @Description("""
        The Seqera Platform compute environment ID. When specified, the scheduler resolves
        the compute environment directly by this ID instead of listing all workspace CEs.
        Used as a fallback when the workflow launch does not include a CE reference.
    """)
    final String computeEnvId

    /* required by config scope -- do not remove */

    ExecutorOpts() {}

    ExecutorOpts(Map opts) {
        this.retryPolicy = new RetryOpts(opts.retryPolicy as Map ?: Map.of())
        this.endpoint = opts.endpoint as String
        if (!endpoint)
            throw new IllegalArgumentException("Missing Seqera endpoint - make sure to specify 'seqera.executor.endpoint' settings")

        this.provider = opts.provider as String
        this.region = opts.region as String
        this.keyPairName = opts.keyPairName as String
        this.batchFlushInterval = opts.batchFlushInterval
            ? Duration.of(opts.batchFlushInterval as String)
            : Duration.of('1 sec')
        // machine requirement settings
        this.machineRequirement = new MachineRequirementOpts(opts.machineRequirement as Map ?: Map.of())
        this.autoLabels = parseAutoLabels(opts.get('autoLabels'))
        // prediction model
        this.predictionModel = opts.predictionModel as String ?: null
        // custom task environment variables
        this.taskEnvironment = opts.taskEnvironment as Map<String, String>
        // compute environment ID
        this.computeEnvId = opts.computeEnvId as String
    }

    RetryOpts retryOpts() {
        this.retryPolicy
    }

    String getEndpoint() {
        return endpoint
    }

    String getProvider() {
        return provider
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

    Set<String> getAutoLabels() {
        return autoLabels
    }

    protected static Set<String> parseAutoLabels(Object value) {
        if( value == null || value == false )
            return Collections.<String>emptySet()
        if( value == true )
            return VALID_AUTO_LABELS
        List<String> raw
        if( value instanceof CharSequence )
            raw = value.toString().tokenize(',').collect { String s -> s.trim() }.findAll { String s -> s }
        else if( value instanceof List )
            raw = ((List) value).collect { it?.toString()?.trim() }.findAll { String s -> s } as List<String>
        else
            throw new IllegalArgumentException("Invalid 'seqera.executor.autoLabels' value '${value}' - expected true, false, a list, or a comma-separated string")
        final invalid = raw.findAll { String s -> !(s in VALID_AUTO_LABELS) }
        if( invalid )
            throw new IllegalArgumentException("Invalid 'seqera.executor.autoLabels' name(s) ${invalid} - valid names are: ${VALID_AUTO_LABELS.join(', ')}")
        return Collections.unmodifiableSet(new LinkedHashSet<>(raw))
    }

    String getPredictionModel() {
        return predictionModel
    }

    Map<String, String> getTaskEnvironment() {
        return taskEnvironment
    }

    String getComputeEnvId() {
        return computeEnvId
    }
}
