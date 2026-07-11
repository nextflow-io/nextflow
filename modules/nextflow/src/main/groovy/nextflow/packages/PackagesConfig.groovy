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

package nextflow.packages

import groovy.transform.CompileStatic
import nextflow.config.schema.ConfigOption
import nextflow.config.schema.ConfigScope
import nextflow.config.schema.ScopeName
import nextflow.script.dsl.Description

/**
 * Model the configuration of the unified `package` directive
 */
@ScopeName("packages")
@Description("""
    The `packages` scope controls the behavior of the `package` process directive
    (requires the `nextflow.preview.package` feature flag).
""")
@CompileStatic
class PackagesConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        The default package provider used when a `package` directive does not specify one (default: `conda`).
    """)
    final String provider

    @ConfigOption
    @Description("""
        Auto-detect a provider manifest file (e.g. `environment.yml`, `requirements.txt`) in the
        process module directory when a process does not declare a `package` directive (default: `true`).
    """)
    final boolean autoDetect

    /* required by extension point -- do not remove */
    PackagesConfig() {}

    PackagesConfig(Map opts) {
        provider = opts.provider != null ? opts.provider as String : 'conda'
        autoDetect = opts.autoDetect != null ? opts.autoDetect as boolean : true
    }
}
