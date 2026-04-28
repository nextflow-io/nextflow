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

@ScopeName("appleContainer")
@Description("""
    The `appleContainer` scope controls how [Apple container](https://github.com/apple/container) runtime is used by Nextflow.
""")
@Slf4j
@CompileStatic
@EqualsAndHashCode
class AppleContainerConfig implements ConfigScope, ContainerConfig {

    @ConfigOption
    @Description("""
        Enable Apple container execution (default: `false`).
    """)
    boolean enabled

    @ConfigOption
    @Description("""
        Specify additional options supported by the `container` CLI i.e. `container [OPTIONS] run`.
    """)
    final String engineOptions

    @ConfigOption
    @Description("""
        Comma separated list of environment variable names to be included in the container environment.
    """)
    final List<String> envWhitelist

    @ConfigOption(types=[String,Boolean])
    @Description("""
    """)
    final Object kill

    @ConfigOption
    @Description("""
        The registry from where container images are pulled. It should NOT include the protocol prefix i.e. `http://`.
    """)
    final String registry

    @ConfigOption
    @Description("""
        Clean up the container after the execution (default: `true`).
    """)
    final boolean remove

    @ConfigOption
    @Description("""
        Specify extra command line options supported by the `container run` command.
    """)
    final String runOptions

    @ConfigOption
    @Description("""
        Mounts a path of your choice as the `/tmp` directory in the container.
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
    AppleContainerConfig() {}

    AppleContainerConfig(Map opts) {
        enabled = opts.enabled as boolean
        engineOptions = opts.engineOptions
        envWhitelist = ContainerHelper.parseEnvWhitelist(opts.envWhitelist)
        kill = opts.kill != null ? opts.kill : true
        registry = opts.registry
        remove = opts.remove != null ? opts.remove as boolean : true
        runOptions = opts.runOptions
        temp = opts.temp
        tty = opts.tty as boolean
        writableInputMounts = opts.writableInputMounts != null ? opts.writableInputMounts as boolean : true
    }

    @Override
    String getEngine() {
        return 'apple-container'
    }

    @Override
    boolean canRunOciImage() {
        return true
    }
}
