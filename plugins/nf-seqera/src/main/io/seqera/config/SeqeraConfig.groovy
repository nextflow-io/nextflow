/*
 * Copyright 2013-2025, Seqera Labs
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

package io.seqera.config

import groovy.transform.CompileStatic
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description
import nextflow.util.Duration

/**
 * Configuration for the Seqera executor.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ScopeName("seqera")
@Description("""
    The `seqera` scope provides configuration for the Seqera compute executor.
""")
@CompileStatic
class SeqeraConfig implements ConfigScope {

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

    @ConfigOption
    @Description("""
        Machine/infrastructure requirements for session tasks.
    """)
    final MachineRequirementOpts machineRequirement

    /* required by config scope -- do not remove */
    SeqeraConfig() {}

    SeqeraConfig(Map opts) {
        this.retryPolicy = new RetryOpts(opts.retryPolicy as Map ?: Map.of())
        this.endpoint = opts.endpoint as String
        if (!endpoint)
            throw new IllegalArgumentException("Missing Seqera endpoint - make sure to specify 'seqera.endpoint' settings")

        this.region = opts.region as String ?: "eu-central-1"
        this.keyPairName = opts.keyPairName as String
        this.batchFlushInterval = opts.batchFlushInterval
            ? Duration.of(opts.batchFlushInterval as String)
            : Duration.of('1 sec')
        // machine requirement settings
        this.machineRequirement = new MachineRequirementOpts(opts.machineRequirement as Map ?: Map.of())
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
}
