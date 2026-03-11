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
package nextflow.container

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description

@ScopeName("sarus")
@Description("""
    The `sarus` scope controls how [Sarus](https://sarus.readthedocs.io) containers are executed by Nextflow.
""")
@CompileStatic
@EqualsAndHashCode
class SarusConfig implements ConfigScope, ContainerConfig {

    @ConfigOption
    @Description("""
        Execute tasks with Sarus containers (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        Comma-separated list of environment variable names to be included in the container environment.
    """)
    final List<String> envWhitelist

    @ConfigOption
    @Description("""
        Specify extra command line options supported by the `sarus run` command.
    """)
    final String runOptions

    @ConfigOption
    @Description("""
        Allocates a pseudo-tty (default: `false`).
    """)
    final String tty

    @ConfigOption
    @Description("""
    """)
    final boolean verbose

    /* required by extension point -- do not remove */
    SarusConfig() {}

    SarusConfig(Map opts) {
        enabled = opts.enabled as boolean
        envWhitelist = ContainerHelper.parseEnvWhitelist(opts.envWhitelist)
        runOptions = opts.runOptions
        tty = opts.tty
        verbose = opts.verbose as boolean
    }

    @Override
    String getEngine() {
        return 'sarus'
    }

}
