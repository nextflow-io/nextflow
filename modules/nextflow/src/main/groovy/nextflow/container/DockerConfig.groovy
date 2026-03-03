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
package nextflow.container

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.util.logging.Slf4j
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description

@ScopeName("docker")
@Description("""
    The `docker` scope controls how [Docker](https://www.docker.com) containers are executed by Nextflow.
""")
@Slf4j
@CompileStatic
@EqualsAndHashCode
class DockerConfig implements ConfigScope, ContainerConfig {

    @ConfigOption
    @Description("""
        Enable Docker execution (default: `false`).
    """)
    boolean enabled

    @ConfigOption
    @Description("""
        Specify additional options supported by the Docker engine i.e. `docker [OPTIONS]`.
    """)
    final String engineOptions

    @ConfigOption
    @Description("""
        Comma separated list of environment variable names to be included in the container environment.
    """)
    final List<String> envWhitelist

    @ConfigOption
    @Description("""
        Fix ownership of files created by the Docker container (default: `false`).
    """)
    final boolean fixOwnership

    @ConfigOption(types=[String,Boolean])
    @Description("""
    """)
    final Object kill

    @ConfigOption
    @Description("""
        Use command line options removed since Docker 1.10.0 (default: `false`).
    """)
    final boolean legacy

    @ConfigOption
    @Description("""
        Add the specified flags to the volume mounts e.g. `'ro,Z'`.
    """)
    final String mountFlags

    @ConfigOption
    @Description("""
        The registry from where Docker images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.
    """)
    final String registry

    @ConfigOption
    @Description("""
        When `true`, forces the override of the registry name in fully qualified container image names with the registry specified by `docker.registry` (default: `false`).
    """)
    final boolean registryOverride

    @ConfigOption
    @Description("""
        Clean up the container after the execution (default: `true`). See the [Docker documentation](https://docs.docker.com/engine/reference/run/#clean-up---rm) for details.
    """)
    final boolean remove

    @ConfigOption
    @Description("""
        Specify extra command line options supported by the `docker run` command. See the [Docker documentation](https://docs.docker.com/engine/reference/run/) for details.
    """)
    final String runOptions

    @ConfigOption
    @Description("""
        Executes Docker run command as `sudo` (default: `false`).
    """)
    final boolean sudo

    @ConfigOption
    @Description("""
        Mounts a path of your choice as the `/tmp` directory in the container. Use the special value `'auto'` to create a temporary directory each time a container is created.
    """)
    final String temp

    @ConfigOption
    @Description("""
        Allocates a pseudo-tty (default: `false`).
    """)
    final boolean tty

    @ConfigOption
    @Description("""
        When `false`, mount input directories as read-only (default: `true`).
    """)
    final boolean writableInputMounts

    /* required by extension point -- do not remove */
    DockerConfig() {}

    DockerConfig(Map opts) {
        enabled = opts.enabled as boolean
        engineOptions = opts.engineOptions
        envWhitelist = ContainerHelper.parseEnvWhitelist(opts.envWhitelist)
        fixOwnership = opts.fixOwnership as boolean
        kill = opts.kill != null ? opts.kill : true
        legacy = opts.legacy as boolean
        mountFlags = opts.mountFlags
        registry = opts.registry
        registryOverride = opts.registryOverride as boolean
        remove = opts.remove != null ? opts.remove as boolean : true
        runOptions = opts.runOptions
        sudo = opts.sudo as boolean
        temp = opts.temp
        tty = opts.tty as boolean
        writableInputMounts = opts.writableInputMounts != null ? opts.writableInputMounts as boolean : true

        if( opts.userEmulation )
            log.warn1("Config setting `docker.userEmulation` is not supported anymore")
    }

    @Override
    String getEngine() {
        return 'docker'
    }

    @Override
    String getFusionOptions() {
        return '--rm --privileged'
    }

}
