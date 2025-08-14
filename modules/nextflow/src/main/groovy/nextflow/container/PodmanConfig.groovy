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
import groovy.transform.EqualsAndHashCode
import nextflow.config.schema.ConfigOption
import nextflow.config.schema.ConfigScope
import nextflow.config.schema.ScopeName
import nextflow.script.dsl.Description

@ScopeName("podman")
@Description("""
    The `podman` scope controls how [Podman](https://podman.io/) containers are executed by Nextflow.
""")
@CompileStatic
@EqualsAndHashCode
class PodmanConfig implements ConfigScope, ContainerConfig {

    @ConfigOption
    @Description("""
        Execute tasks with Podman containers (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        Specify additional options supported by the Podman engine i.e. `podman [OPTIONS]`.
    """)
    final String engineOptions

    @ConfigOption
    @Description("""
        Comma separated list of environment variable names to be included in the container environment.
    """)
    final List<String> envWhitelist

    @ConfigOption
    @Description("""
    """)
    final Object kill

    @ConfigOption
    @Description("""
        Add the specified flags to the volume mounts e.g. `'ro,Z'`.
    """)
    final String mountFlags

    @ConfigOption
    @Description("""
        The registry from where container images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.
    """)
    final String registry

    @ConfigOption
    @Description("""
        Clean-up the container after the execution (default: `true`).
    """)
    final boolean remove

    @ConfigOption
    @Description("""
        Specify extra command line options supported by the `podman run` command.
    """)
    final String runOptions

    @ConfigOption
    @Description("""
        Mounts a path of your choice as the `/tmp` directory in the container. Use the special value `'auto'` to create a temporary directory each time a container is created.
    """)
    final String temp

    /* required by extension point -- do not remove */
    PodmanConfig() {}

    PodmanConfig(Map opts) {
        enabled = opts.enabled as boolean
        engineOptions = opts.engineOptions
        envWhitelist = ContainerHelper.parseEnvWhitelist(opts.envWhitelist)
        kill = opts.kill != null ? opts.kill : true
        mountFlags = opts.mountFlags
        registry = opts.registry
        remove = opts.remove != null ? opts.remove as boolean : true
        runOptions = opts.runOptions
        temp = opts.temp
    }

    @Override
    String getEngine() {
        return 'podman'
    }

    @Override
    String getFusionOptions() {
        return '--rm --privileged'
    }

}
