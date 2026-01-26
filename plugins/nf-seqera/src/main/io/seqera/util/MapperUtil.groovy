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
import io.seqera.sched.api.schema.v1a1.MachineRequirement
import io.seqera.sched.api.schema.v1a1.PriceModel as SchedPriceModel
import io.seqera.sched.api.schema.v1a1.ProvisioningModel
import nextflow.cloud.types.PriceModel

/**
 * Utility class to map Nextflow config objects to Sched API schema objects.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class MapperUtil {

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
        final arch = taskArch ?: opts?.arch
        final provisioning = opts?.provisioning
        final maxSpotAttempts = opts?.maxSpotAttempts
        final machineFamilies = opts?.machineFamilies
        // return null if no settings
        if (!arch && !provisioning && !maxSpotAttempts && !machineFamilies)
            return null
        new MachineRequirement()
            .arch(arch)
            .provisioning(toProvisioningModel(provisioning))
            .maxSpotAttempts(maxSpotAttempts)
            .machineFamilies(machineFamilies)
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
