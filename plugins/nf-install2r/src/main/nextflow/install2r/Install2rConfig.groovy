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

package nextflow.install2r

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.config.schema.ConfigOption
import nextflow.config.schema.ConfigScope
import nextflow.config.schema.ScopeName
import nextflow.script.dsl.Description
import nextflow.util.Duration

/**
 * Model install2.r configuration
 */
@ScopeName("install2r")
@Description("""
    The `install2r` scope controls the installation of R/CRAN packages by the install2.r (littler) package manager.
""")
@CompileStatic
class Install2rConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        Execute tasks with install2.r R libraries (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        The path where install2.r R libraries are stored. It should be accessible from all compute nodes when using a shared file system.
    """)
    final String cacheDir

    @ConfigOption
    @Description("""
        Extra command line options for the `install2.r` command. See the littler documentation for more information.
    """)
    final String installOptions

    @ConfigOption
    @Description("""
        The amount of time to wait for the R library to be created before failing (default: `20 min`).
    """)
    final Duration createTimeout

    @ConfigOption
    @Description("""
        The CRAN repository used by `install2.r` to resolve packages (default: `https://cloud.r-project.org`).
    """)
    final String repos

    /* required by extension point -- do not remove */
    Install2rConfig() {}

    Install2rConfig(Map opts, Map<String, String> env) {
        enabled = opts.enabled != null
            ? opts.enabled as boolean
            : (env.NXF_INSTALL2R_ENABLED?.toString() == 'true')
        cacheDir = opts.cacheDir
        installOptions = opts.installOptions
        createTimeout = opts.createTimeout as Duration ?: Duration.of('20min')
        repos = opts.repos ?: 'https://cloud.r-project.org'
    }

    boolean isEnabled() {
        enabled
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

    String repos() {
        repos
    }
}
