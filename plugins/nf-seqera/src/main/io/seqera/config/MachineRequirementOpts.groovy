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
import nextflow.script.dsl.Description

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
        List of acceptable EC2 instance families (e.g., `['m5', 'c5', 'r5']`).
    """)
    final List<String> machineFamilies

    /* required by config scope -- do not remove */
    MachineRequirementOpts() {}

    MachineRequirementOpts(Map opts) {
        this.arch = opts.arch as String
        this.provisioning = opts.provisioning as String
        this.maxSpotAttempts = opts.maxSpotAttempts as Integer
        this.machineFamilies = opts.machineFamilies as List<String>
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

    List<String> getMachineFamilies() {
        return machineFamilies
    }
}