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
import nextflow.platform.AutoLabels
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

    @Description("""
        HTTP client settings for requests to the Seqera scheduler service.
    """)
    final HttpClientOpts httpClient

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
        Execution strategy within the chosen provider.
        Narrows compute-environment selection when a provider offers multiple strategies
        (e.g. AWS supports `ecs` and `vm`). When omitted, the provider's canonical default
        strategy is used (AWS → `ecs`).
    """)
    final String strategy

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
        The interval for batching task submissions (default: `5 sec`).
    """)
    final Duration batchFlushInterval

    @Description("""
        Machine/infrastructure requirements for session tasks.
    """)
    final MachineRequirementOpts machineRequirement

    @Deprecated
    @ConfigOption
    @Description("""
        DEPRECATED: use `tower.autoLabels` instead. This option is honoured for backward
        compatibility and takes precedence over `tower.autoLabels` when specified, including
        when set to `false`.

        Automatically attach workflow metadata labels (with the `nextflow.io/` and
        `seqera.io/platform/` prefixes) to the compute resources. Accepts:
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

    @Description("""
        Scheduling requirements applied to this run by the Seqera scheduler.
    """)
    final SchedulingRequirementOpts schedulingRequirement

    @ConfigOption
    @Description("""
        Custom environment variables to apply to all tasks submitted by the Seqera executor.
        These are merged with the Fusion environment variables, with Fusion variables taking precedence.
    """)
    final Map<String, String> taskEnvironment

    @ConfigOption
    @Description("""
        Backend-specific provider configuration merged into the compute cluster's backend
        properties (for cluster isolation). When omitted, the backend falls back to its
        environment variable configuration.
    """)
    final Map<String, String> providerConfig

    @ConfigOption
    @Description("""
        The Seqera Platform compute environment ID. When specified, the scheduler resolves
        the compute environment directly by this ID instead of listing all workspace CEs.
        Used as a fallback when the workflow launch does not include a CE reference.
    """)
    final String computeEnvId

    @ConfigOption
    @Description("""
        Enable on-demand interactive shell access (e.g. SSH) to this run's task containers
        (VM and local backends). When `true`, a running task can be reached with
        `sched task ssh <task-id>` (or a plain `ssh <task-id>@<scheduler>`), and the
        connection survives task completion. Default: `false`.
    """)
    final boolean shellEnabled

    /* required by config scope -- do not remove */

    ExecutorOpts() {}

    ExecutorOpts(Map opts) {
        this.retryPolicy = new RetryOpts(opts.retryPolicy as Map ?: Map.of())
        this.httpClient = new HttpClientOpts(opts.httpClient as Map ?: Map.of())
        this.endpoint = opts.endpoint as String
        if (!endpoint)
            throw new IllegalArgumentException("Missing Seqera endpoint - make sure to specify 'seqera.executor.endpoint' settings")

        this.provider = opts.provider as String
        this.strategy = opts.strategy as String
        this.region = opts.region as String
        this.keyPairName = opts.keyPairName as String
        this.batchFlushInterval = opts.batchFlushInterval
            ? Duration.of(opts.batchFlushInterval as String)
            : Duration.of('5 sec')
        // machine requirement settings
        this.machineRequirement = new MachineRequirementOpts(opts.machineRequirement as Map ?: Map.of())
        // note the labels are resolved session-wide by `Session#getAutoResourceLabels` -- parsing
        // them here still validates the setting and rejects an unknown name at config load time
        this.autoLabels = AutoLabels.parse(opts.get('autoLabels'), 'seqera.executor.autoLabels')
        // prediction model
        this.predictionModel = opts.predictionModel as String ?: null
        // scheduling requirements (e.g. per-user vCPU cap)
        this.schedulingRequirement = new SchedulingRequirementOpts(opts.schedulingRequirement as Map ?: Map.of())
        // custom task environment variables
        this.taskEnvironment = opts.taskEnvironment as Map<String, String>
        // backend-specific provider configuration
        this.providerConfig = opts.providerConfig as Map<String, String>
        // compute environment ID
        this.computeEnvId = opts.computeEnvId as String
        // on-demand shell access to task containers (default false)
        this.shellEnabled = opts.shellEnabled as boolean
    }

    RetryOpts retryOpts() {
        this.retryPolicy
    }

    HttpClientOpts httpOpts() {
        this.httpClient
    }

    String getEndpoint() {
        return endpoint
    }

    String getProvider() {
        return provider
    }

    String getStrategy() {
        return strategy
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

    String getPredictionModel() {
        return predictionModel
    }

    SchedulingRequirementOpts getSchedulingRequirement() {
        return schedulingRequirement
    }

    Map<String, String> getTaskEnvironment() {
        return taskEnvironment
    }

    Map<String, String> getProviderConfig() {
        return providerConfig
    }

    String getComputeEnvId() {
        return computeEnvId
    }

    boolean getShellEnabled() {
        return shellEnabled
    }
}
