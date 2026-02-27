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
import nextflow.util.MemoryUnit

/**
 * Machine/infrastructure requirements configuration options.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class MachineRequirementOpts implements ConfigScope {

    @ConfigOption
    @Description("""
        The CPU architecture for task execution (e.g., `x86_64`, `arm64`).
    """)
    final String arch

    @ConfigOption
    @Description("""
        The instance provisioning mode: `spot`, `ondemand`, or `spotFirst`.
    """)
    final String provisioning

    @ConfigOption
    @Description("""
        Maximum number of spot retry attempts before falling back to on-demand.
        Only used when provisioning is `spot` or `spotFirst`.
    """)
    final Integer maxSpotAttempts

    @ConfigOption
    @Description("""
        List of acceptable machine type patterns. Supports exact types (e.g., `t3.small`),
        family prefixes (e.g., `m5` matches all m5 sizes), and glob wildcards (e.g., `t*.small`).
    """)
    final List<String> machineTypes

    @ConfigOption
    @Description("""
        The EBS volume type for task scratch disk (e.g., `ebs/gp3`, `ebs/io1`).
        Default: `ebs/gp3`.
    """)
    final String diskType

    @ConfigOption
    @Description("""
        The throughput in MiB/s for gp3 volumes (125-1000).
        Default: 325 (Fusion recommended).
    """)
    final Integer diskThroughputMiBps

    @ConfigOption
    @Description("""
        The IOPS for io1/io2/gp3 volumes. Required for io1/io2.
    """)
    final Integer diskIops

    @ConfigOption
    @Description("""
        Enable KMS encryption for the EBS volume.
        Default: false.
    """)
    final Boolean diskEncrypted

    @ConfigOption
    @Description("""
        The disk allocation strategy: `task` or `node`.
        - `task`: Per-task EBS volume created at task launch (default)
        - `node`: Per-node instance storage attached at cluster level
    """)
    final String diskAllocation

    @ConfigOption
    @Description("""
        The disk size for session-level storage (e.g., `100.GB`).
    """)
    final MemoryUnit diskSize

    @ConfigOption
    @Description("""
        The ECS capacity provider mode: `managed` (default) or `asg`.
        - `managed`: ECS Managed Instances
        - `asg`: Auto Scaling Group-backed capacity provider
    """)
    final String capacityMode

    /* required by config scope -- do not remove */
    MachineRequirementOpts() {}

    MachineRequirementOpts(Map opts) {
        this.arch = opts.arch as String
        this.provisioning = opts.provisioning as String
        this.maxSpotAttempts = opts.maxSpotAttempts as Integer
        this.machineTypes = (opts.machineTypes ?: opts.machineFamilies) as List<String>
        this.diskType = opts.diskType as String
        this.diskThroughputMiBps = opts.diskThroughputMiBps as Integer
        this.diskIops = opts.diskIops as Integer
        this.diskEncrypted = opts.diskEncrypted as Boolean
        this.diskAllocation = opts.diskAllocation as String
        this.diskSize = opts.diskSize instanceof MemoryUnit
            ? opts.diskSize as MemoryUnit
            : (opts.diskSize ? MemoryUnit.of(opts.diskSize as String) : null)
        this.capacityMode = opts.capacityMode as String
    }

    String getArch() {
        return arch
    }

    String getProvisioning() {
        return provisioning
    }

    Integer getMaxSpotAttempts() {
        return maxSpotAttempts
    }

    List<String> getMachineTypes() {
        return machineTypes
    }

    String getDiskType() {
        return diskType
    }

    Integer getDiskThroughputMiBps() {
        return diskThroughputMiBps
    }

    Integer getDiskIops() {
        return diskIops
    }

    Boolean getDiskEncrypted() {
        return diskEncrypted
    }

    String getDiskAllocation() {
        return diskAllocation
    }

    MemoryUnit getDiskSize() {
        return diskSize
    }

    String getCapacityMode() {
        return capacityMode
    }
}
