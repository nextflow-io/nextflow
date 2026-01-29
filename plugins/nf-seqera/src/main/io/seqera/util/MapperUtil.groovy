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

package io.seqera.util

import groovy.transform.CompileStatic
import io.seqera.config.MachineRequirementOpts
import io.seqera.sched.api.schema.v1a1.DiskRequirement
import io.seqera.sched.api.schema.v1a1.MachineRequirement
import io.seqera.sched.api.schema.v1a1.PriceModel as SchedPriceModel
import io.seqera.sched.api.schema.v1a1.ProvisioningModel
import nextflow.cloud.types.PriceModel
import nextflow.util.MemoryUnit

/**
 * Utility class to map Nextflow config objects to Sched API schema objects.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class MapperUtil {

    /** Default EBS volume type - gp3 provides good balance of price and performance */
    static final String DEFAULT_DISK_TYPE = 'ebs/gp3'

    /** Default throughput in MiB/s - Fusion recommended setting for optimal I/O */
    static final int DEFAULT_DISK_THROUGHPUT_MIBPS = 325

    /** Supported EBS volume types */
    static final Set<String> SUPPORTED_DISK_TYPES = Set.of(
        'ebs/gp3',   // General purpose SSD (default)
        'ebs/gp2',   // General purpose SSD (legacy)
        'ebs/io1',   // Provisioned IOPS SSD
        'ebs/io2',   // Provisioned IOPS SSD (higher durability)
        'ebs/st1',   // Throughput optimized HDD
        'ebs/sc1'    // Cold HDD
    )

    /**
     * Maps MachineRequirementOpts to MachineRequirement API object.
     *
     * @param opts the config options (can be null)
     * @return the MachineRequirement API object, or null if opts is null or has no settings
     */
    static MachineRequirement toMachineRequirement(MachineRequirementOpts opts) {
        if (!opts)
            return null
        if (!opts.arch && !opts.provisioning && !opts.maxSpotAttempts && !opts.machineFamilies)
            return null
        new MachineRequirement()
            .arch(opts.arch)
            .provisioning(toProvisioningModel(opts.provisioning))
            .maxSpotAttempts(opts.maxSpotAttempts)
            .machineFamilies(opts.machineFamilies)
    }

    /**
     * Maps MachineRequirementOpts to MachineRequirement API object, merging with task arch.
     * Task arch overrides config arch if specified.
     *
     * @param opts the config options (can be null)
     * @param taskArch the task container platform/arch (can be null)
     * @return the MachineRequirement API object, or null if no settings
     */
    static MachineRequirement toMachineRequirement(MachineRequirementOpts opts, String taskArch) {
        return toMachineRequirement(opts, taskArch, null)
    }

    /**
     * Maps MachineRequirementOpts to MachineRequirement API object, merging with task arch and disk.
     * Task arch overrides config arch if specified.
     *
     * @param opts the config options (can be null)
     * @param taskArch the task container platform/arch (can be null)
     * @param diskSize the disk size from task config (can be null)
     * @return the MachineRequirement API object, or null if no settings
     */
    static MachineRequirement toMachineRequirement(MachineRequirementOpts opts, String taskArch, MemoryUnit diskSize) {
        final arch = taskArch ?: opts?.arch
        final provisioning = opts?.provisioning
        final maxSpotAttempts = opts?.maxSpotAttempts
        final machineFamilies = opts?.machineFamilies
        final diskReq = toDiskRequirement(diskSize, opts)
        // return null if no settings
        if (!arch && !provisioning && !maxSpotAttempts && !machineFamilies && !diskReq)
            return null
        new MachineRequirement()
            .arch(arch)
            .provisioning(toProvisioningModel(provisioning))
            .maxSpotAttempts(maxSpotAttempts)
            .machineFamilies(machineFamilies)
            .disk(diskReq)
    }

    /**
     * Maps a disk size to DiskRequirement API object.
     * Uses config options if provided, otherwise defaults to Fusion recommended settings:
     * EBS gp3 volume with 325 MiB/s throughput.
     *
     * @param diskSize the disk size (can be null)
     * @param opts the machine requirement options with disk settings (can be null)
     * @return the DiskRequirement API object, or null if diskSize is null or zero
     */
    static DiskRequirement toDiskRequirement(MemoryUnit diskSize, MachineRequirementOpts opts=null) {
        if (!diskSize || diskSize.toGiga() <= 0)
            return null
        // Use config values or Fusion recommended defaults
        final type = opts?.diskType ?: DEFAULT_DISK_TYPE
        // Validate disk type is supported
        if (!SUPPORTED_DISK_TYPES.contains(type)) {
            throw new IllegalArgumentException("Invalid disk type: ${type}. Supported types: ${SUPPORTED_DISK_TYPES.join(', ')}")
        }
        final throughput = opts?.diskThroughputMiBps ?: DEFAULT_DISK_THROUGHPUT_MIBPS
        final iops = opts?.diskIops
        final encrypted = opts?.diskEncrypted ?: false

        def req = new DiskRequirement()
            .sizeGiB(diskSize.toGiga() as Integer)
            .type(type)
            .encrypted(encrypted)
        // Only set throughput for gp3 volumes
        if (type == DEFAULT_DISK_TYPE) {
            req.throughputMiBps(throughput)
        }
        // Set IOPS if provided
        if (iops) {
            req.iops(iops)
        }
        return req
    }

    /**
     * Maps a provisioning string to ProvisioningModel enum.
     *
     * @param value the provisioning string (spot, ondemand, spotFirst)
     * @return the ProvisioningModel enum value, or null if value is null
     */
    static ProvisioningModel toProvisioningModel(String value) {
        value ? ProvisioningModel.fromValue(value) : null
    }

    /**
     * Maps Sched API PriceModel to Nextflow PriceModel.
     *
     * @param schedPriceModel the Sched API price model
     * @return the Nextflow PriceModel, or null if input is null or unknown
     */
    static PriceModel toPriceModel(SchedPriceModel schedPriceModel) {
        if (schedPriceModel == null)
            return null
        switch (schedPriceModel) {
            case SchedPriceModel.SPOT:
                return PriceModel.spot
            case SchedPriceModel.STANDARD:
                return PriceModel.standard
            default:
                return null
        }
    }

}
