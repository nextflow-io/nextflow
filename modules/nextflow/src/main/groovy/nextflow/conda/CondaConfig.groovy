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

import java.nio.file.Path

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
""")
@CompileStatic
class CondaConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        Execute tasks with Conda environments (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        The path where Conda environments are stored. It should be accessible from all compute nodes when using a shared file system.
    """)
    final String cacheDir

    @ConfigOption
    @Description("""
        The list of Conda channels that can be used to resolve Conda packages. Channel priority decreases from left to right.
    """)
    final List<String> channels

    @ConfigOption
    @Description("""
        Extra command line options for the `conda create` command. See the [Conda documentation](https://docs.conda.io/projects/conda/en/latest/commands/create.html) for more information.
    """)
    final String createOptions

    @ConfigOption
    @Description("""
        The amount of time to wait for the Conda environment to be created before failing (default: `20 min`).
    """)
    final Duration createTimeout

    @ConfigOption
    @Description("""
        Use [Mamba](https://github.com/mamba-org/mamba) instead of `conda` to create Conda environments (default: `false`).
    """)
    final boolean useMamba

    @ConfigOption
    @Description("""
        Use [Micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) instead of `conda` to create Conda environments (default: `false`).
    """)
    final boolean useMicromamba

    /* required by extension point -- do not remove */
    CondaConfig() {}

    CondaConfig(Map opts, Map<String, String> env) {
        enabled = opts.enabled != null
            ? opts.enabled as boolean
            : (env.NXF_CONDA_ENABLED?.toString() == 'true')
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
        if( value instanceof CharSequence )
            return value.tokenize(',').collect(it -> it.trim())
        throw new IllegalArgumentException("Unexpected conda.channels value: $value")
    }

    Duration createTimeout() {
        createTimeout
    }

    String createOptions() {
        createOptions
    }

    Path cacheDir() {
        cacheDir as Path
    }

    boolean useMamba() {
        useMamba
    }

    boolean useMicromamba() {
        useMicromamba
    }
}
