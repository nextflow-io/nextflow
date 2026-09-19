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

package nextflow.pixi

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description
import nextflow.util.Duration

/**
 * Model Pixi configuration
 *
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
@ScopeName("pixi")
@Description("""
    The `pixi` scope controls the creation of package environments by the Pixi package manager.
""")
@CompileStatic
class PixiConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        The path where Pixi environments are stored. It should be accessible from all compute nodes when using a shared file system.
    """)
    final String cacheDir

    @ConfigOption
    @Description("""
        The list of conda channels used to resolve packages given as a package list (default: `['conda-forge']`).
    """)
    final List<String> channels

    @ConfigOption
    @Description("""
        Extra command line options for the `pixi` commands used to create environments. See the [Pixi documentation](https://pixi.sh/latest/) for more information.
    """)
    final String createOptions

    @ConfigOption
    @Description("""
        The amount of time to wait for the Pixi environment to be created before failing (default: `20 min`).
    """)
    final Duration createTimeout

    /* required by extension point -- do not remove */
    PixiConfig() {}

    PixiConfig(Map opts, Map<String, String> env) {
        cacheDir = opts.cacheDir
        channels = opts.channels != null ? (opts.channels as List).collect { it.toString() } : ['conda-forge']
        createOptions = opts.createOptions
        createTimeout = opts.createTimeout as Duration ?: Duration.of('20min')
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

    List<String> channels() { channels }
}
