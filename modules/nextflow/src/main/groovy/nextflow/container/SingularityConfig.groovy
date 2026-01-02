/*
 * Copyright 2024-2025, Seqera Labs
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
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description
import nextflow.util.Duration

@ScopeName("singularity")
@Description("""
    The `singularity` scope controls how [Singularity](https://sylabs.io/singularity/) containers are executed by Nextflow.
""")
@CompileStatic
class SingularityConfig implements ConfigScope, ContainerConfig {

    @ConfigOption
    @Description("""
        Automatically mount host paths in the executed container (default: `true`). It requires the `user bind control` feature to be enabled in your Singularity installation.
    """)
    final Boolean autoMounts

    @ConfigOption
    @Description("""
        The directory where remote Singularity images are stored. When using a compute cluster, it must be a shared folder accessible to all compute nodes.
    """)
    final String cacheDir

    @ConfigOption
    @Description("""
        Execute tasks with Singularity containers (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        Specify additional options supported by the Singularity engine i.e. `singularity [OPTIONS]`.
    """)
    final String engineOptions

    @ConfigOption
    @Description("""
        Comma separated list of environment variable names to be included in the container environment.
    """)
    final List<String> envWhitelist

    @ConfigOption
    @Description("""
        Directory where remote Singularity images are retrieved. When using a computing cluster it must be a shared folder accessible to all compute nodes.
    """)
    final String libraryDir

    @ConfigOption
    @Description("""
        Pull the Singularity image with http protocol (default: `false`).
    """)
    final boolean noHttps

    @ConfigOption
    @Description("""
        When enabled, OCI (and Docker) container images are pull and converted to a SIF image file format implicitly by the Singularity run command, instead of Nextflow (default: `false`).
    """)
    final boolean ociAutoPull

    @ConfigOption
    @Description("""
        Enable OCI-mode, that allows running native OCI compliant container image with Singularity using `crun` or `runc` as low-level runtime (default: `false`).
    """)
    final Boolean ociMode

    @ConfigOption
    @Description("""
        The amount of time the Singularity pull can last, after which the process is terminated (default: `20 min`).
    """)
    final Duration pullTimeout

    @ConfigOption
    @Description("""
        The registry from where Docker images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.
    """)
    final String registry

    @ConfigOption
    @Description("""
        Specify extra command line options supported by `singularity exec`.
    """)
    final String runOptions

    /* required by extension point -- do not remove */
    SingularityConfig() {}

    SingularityConfig(Map opts) {
        autoMounts = opts.autoMounts as Boolean
        cacheDir = opts.cacheDir
        enabled = opts.enabled as boolean
        engineOptions = opts.engineOptions
        envWhitelist = ContainerHelper.parseEnvWhitelist(opts.envWhitelist)
        libraryDir = opts.libraryDir
        noHttps = opts.noHttps as boolean
        ociAutoPull = opts.ociAutoPull as boolean
        ociMode = opts.ociMode as Boolean
        pullTimeout = opts.pullTimeout as Duration ?: Duration.of('20min')
        registry = opts.registry
        runOptions = opts.runOptions
    }

    @Override
    String getEngine() {
        return 'singularity'
    }

    @Override
    boolean canRunOciImage() {
        return ociMode || ociAutoPull
    }

    @Override
    String getFusionOptions() {
        return ociMode ? '-B /dev/fuse' : null
    }

}
