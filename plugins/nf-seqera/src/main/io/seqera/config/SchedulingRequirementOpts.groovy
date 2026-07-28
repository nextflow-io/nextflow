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

/**
 * Scheduling requirement configuration options.
 *
 * @author Jon Martí <jonathan.marti@seqera.io>
 */
@CompileStatic
class SchedulingRequirementOpts implements ConfigScope {

    @ConfigOption
    @Description("""
        Maximum concurrent vCPUs a single user may consume across all of that user's runs on this
        compute environment. Applied independently to each user sharing the CE. When not set, no
        CPU limit is enforced. A task whose own request exceeds the cap fails outright, so keep the
        cap at least as large as the largest per-task `cpus` request (or clamp it with
        `process.resourceLimits.cpus`).
    """)
    final Integer maxCpusPerUser

    /* required by config scope -- do not remove */
    SchedulingRequirementOpts() {}

    SchedulingRequirementOpts(Map opts) {
        this.maxCpusPerUser = opts.maxCpusPerUser as Integer
        if( maxCpusPerUser != null && maxCpusPerUser <= 0 )
            throw new IllegalArgumentException("Invalid 'seqera.executor.schedulingRequirement.maxCpusPerUser' value '${maxCpusPerUser}' - expected a positive integer")
    }

    Integer getMaxCpusPerUser() {
        return maxCpusPerUser
    }
}
