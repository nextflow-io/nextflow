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

package nextflow.uv

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description
import nextflow.util.Duration

/**
 * Model uv configuration
 *
 * @author Evan Floden
 */
@ScopeName("uv")
@Description("""
    The `uv` scope controls the creation of Python virtual environments by the uv package manager.
""")
@CompileStatic
class UvConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        Execute tasks with uv virtual environments (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        The path where uv virtual environments are stored. It should be accessible from all compute nodes when using a shared file system.
    """)
    final String cacheDir

    @ConfigOption
    @Description("""
        Extra command line options for the `uv pip install` command. See the [uv documentation](https://docs.astral.sh/uv/) for more information.
    """)
    final String installOptions

    @ConfigOption
    @Description("""
        The amount of time to wait for the uv environment to be created before failing (default: `20 min`).
    """)
    final Duration createTimeout

    @ConfigOption
    @Description("""
        The Python version to use when creating virtual environments (e.g. `3.12`). If not specified, uv will use its default Python resolution.
    """)
    final String pythonVersion

    /* required by extension point -- do not remove */
    UvConfig() {}

    UvConfig(Map opts, Map<String, String> env) {
        enabled = opts.enabled != null
            ? opts.enabled as boolean
            : (env.NXF_UV_ENABLED?.toString() == 'true')
        cacheDir = opts.cacheDir
        installOptions = opts.installOptions
        createTimeout = opts.createTimeout as Duration ?: Duration.of('20min')
        pythonVersion = opts.pythonVersion
    }

    Duration createTimeout() {
        createTimeout
    }

    String installOptions() {
        installOptions
    }

    Path cacheDir() {
        cacheDir as Path
    }

    String pythonVersion() {
        pythonVersion
    }
}
