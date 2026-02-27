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

@ScopeName("charliecloud")
@Description("""
    The `charliecloud` scope controls how [Charliecloud](https://hpc.github.io/charliecloud/) containers are executed by Nextflow.
""")
@CompileStatic
class CharliecloudConfig implements ConfigScope, ContainerConfig {

    @ConfigOption
    @Description("""
        The directory where remote Charliecloud images are stored. When using a computing cluster it must be a shared folder accessible to all compute nodes.
    """)
    final String cacheDir

    @ConfigOption
    @Description("""
        Execute tasks with Charliecloud containers (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        Comma separated list of environment variable names to be included in the container environment.
    """)
    final List<String> envWhitelist

    @ConfigOption
    @Description("""
        The amount of time the Charliecloud pull can last, exceeding which the process is terminated (default: `20 min`).
    """)
    final Duration pullTimeout

    @ConfigOption
    @Description("""
        The registry from where images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.
    """)
    final String registry

    @ConfigOption
    @Description("""
        Specify extra command line options supported by the `ch-run` command.
    """)
    final String runOptions

    @ConfigOption
    @Description("""
        Mounts a path of your choice as the `/tmp` directory in the container. Use the special value `'auto'` to create a temporary directory each time a container is created.
    """)
    final String temp

    @ConfigOption
    @Description("""
        When `false`, mount input directories as read-only (default: `true`).
    """)
    final boolean writableInputMounts

    @ConfigOption
    @Description("""
        Run containers from storage in writeable mode using overlayfs (default: `true`).
    """)
    final boolean writeFake

    /* required by extension point -- do not remove */
    CharliecloudConfig() {}

    CharliecloudConfig(Map opts) {
        cacheDir = opts.cacheDir
        enabled = opts.enabled as boolean
        envWhitelist = ContainerHelper.parseEnvWhitelist(opts.envWhitelist)
        pullTimeout = opts.pullTimeout as Duration ?: Duration.of('20min')
        registry = opts.registry
        runOptions = opts.runOptions
        temp = opts.temp
        writableInputMounts = opts.writableInputMounts != null ? opts.writableInputMounts as boolean : true
        writeFake = opts.writeFake != null ? opts.writeFake as boolean : true
    }

    @Override
    String getEngine() {
        return 'charliecloud'
    }

}
