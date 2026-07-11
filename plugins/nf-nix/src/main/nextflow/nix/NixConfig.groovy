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

package nextflow.nix

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.config.schema.ConfigOption
import nextflow.config.schema.ConfigScope
import nextflow.config.schema.ScopeName
import nextflow.script.dsl.Description
import nextflow.util.Duration

/**
 * Model Nix configuration
 */
@ScopeName("nix")
@Description("""
    The `nix` scope controls the creation of environments by the Nix package manager.
""")
@CompileStatic
class NixConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        Execute tasks with Nix environments (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        The path where Nix environments are stored. It should be accessible from all compute nodes when using a shared file system.
    """)
    final String cacheDir

    @ConfigOption
    @Description("""
        Extra command line options appended to the Nix install command. See the [Nix documentation](https://nixos.org/manual/nix/stable/) for more information.
    """)
    final String installOptions

    @ConfigOption
    @Description("""
        The amount of time to wait for the Nix environment to be created before failing (default: `20 min`).
    """)
    final Duration createTimeout

    @ConfigOption
    @Description("""
        The flake reference used to resolve bare package names (default: `nixpkgs`).
    """)
    final String flakeRef

    /* required by extension point -- do not remove */
    NixConfig() {}

    NixConfig(Map opts, Map<String, String> env) {
        enabled = opts.enabled != null
            ? opts.enabled as boolean
            : (env.NXF_NIX_ENABLED?.toString() == 'true')
        cacheDir = opts.cacheDir
        installOptions = opts.installOptions
        createTimeout = opts.createTimeout as Duration ?: Duration.of('20min')
        flakeRef = opts.flakeRef ?: 'nixpkgs'
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

    String flakeRef() {
        flakeRef
    }
}
