/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.conda

import groovy.transform.CompileStatic
import nextflow.config.schema.ConfigOption
import nextflow.config.schema.ConfigScope
import nextflow.config.schema.ScopeName
import nextflow.script.dsl.Description
import nextflow.util.Duration

/**
 * Model Conda configuration
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ScopeName("conda")
@Description("""
    The `conda` scope controls the creation of Conda environments by the Conda package manager.

    [Read more](https://nextflow.io/docs/latest/reference/config.html#conda)
""")
@CompileStatic
class CondaConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        Enable Conda execution (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        The path where Conda environments are stored.
    """)
    final String cacheDir

    @ConfigOption
    @Description("""
        The Conda channels that can be used to resolve Conda packages.
    """)
    final List<String> channels

    @ConfigOption
    @Description("""
        Extra command line options to append to the `conda create` command.
    """)
    final String createOptions

    @ConfigOption
    @Description("""
        The amount of time to wait for the Conda environment to be created before failing (default: `20 min`).
    """)
    final Duration createTimeout

    @ConfigOption
    @Description("""
        When `true`, use `mamba` instead of `conda` to create the Conda environments.

        [Read more](https://github.com/mamba-org/mamba)
    """)
    final boolean useMamba

    @ConfigOption
    @Description("""
        When `true`, use `micromamba` instead of `conda` to create the Conda environments.

        [Read more](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
    """)
    final boolean useMicromamba

    /* required by extension point -- do not remove */
    CondaConfig() {}

    CondaConfig(Map opts, Map<String, String> env) {
        enabled = opts.enabled != null
            ? opts.enabled as boolean
            : (env.NXF_CONDA_ENABLED == 'true')
        cacheDir = opts.cacheDir
        channels = parseChannels(opts.channels)
        createOptions = opts.createOptions
        createTimeout = opts.createTimeout as Duration ?: Duration.of('20min')
        useMamba = opts.useMamba as boolean
        useMicromamba = opts.useMicromamba as boolean

        if( useMamba && useMicromamba )
            throw new IllegalArgumentException("Both conda.useMamba and conda.useMicromamba were enabled -- Please choose only one")
    }

    private List<String> parseChannels(Object value) {
        if( !value )
            return Collections.emptyList()
        if( value instanceof List )
            return value
        throw new IllegalArgumentException("Unexpected conda.channels value: $value")
    }
}
