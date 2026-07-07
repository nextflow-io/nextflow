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

package nextflow.guix

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.config.schema.ConfigOption
import nextflow.config.schema.ConfigScope
import nextflow.config.schema.ScopeName
import nextflow.script.dsl.Description
import nextflow.util.Duration

/**
 * Model GNU Guix configuration
 */
@ScopeName("guix")
@Description("""
    The `guix` scope controls the creation of GNU Guix environments by the Guix package manager.
""")
@CompileStatic
class GuixConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        Execute tasks with GNU Guix environments (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        The path where GNU Guix environments are stored. It should be accessible from all compute nodes when using a shared file system.
    """)
    final String cacheDir

    @ConfigOption
    @Description("""
        Extra command line options for the `guix package --install` command. See the [GNU Guix documentation](https://guix.gnu.org/) for more information.
    """)
    final String installOptions

    @ConfigOption
    @Description("""
        The amount of time to wait for the GNU Guix environment to be created before failing (default: `20 min`).
    """)
    final Duration createTimeout

    /* required by extension point -- do not remove */
    GuixConfig() {}

    GuixConfig(Map opts, Map<String, String> env) {
        enabled = opts.enabled != null
            ? opts.enabled as boolean
            : (env.NXF_GUIX_ENABLED?.toString() == 'true')
        cacheDir = opts.cacheDir
        installOptions = opts.installOptions
        createTimeout = opts.createTimeout as Duration ?: Duration.of('20min')
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
}
