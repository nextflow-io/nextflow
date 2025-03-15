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
package nextflow.config.scopes;

import nextflow.config.schema.ConfigOption;
import nextflow.config.schema.ConfigScope;
import nextflow.script.dsl.Description;

public class DockerConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        Enable Docker execution (default: `false`).
    """)
    public boolean enabled;

    @ConfigOption
    @Description("""
        This attribute can be used to provide any option supported by the Docker engine i.e. `docker [OPTIONS]`.
    """)
    public String engineOptions;

    @ConfigOption
    @Description("""
        Comma separated list of environment variable names to be included in the container environment.
    """)
    public String envWhitelist;

    @ConfigOption
    @Description("""
        Fix ownership of files created by the Docker container.
    """)
    public boolean fixOwnership;

    @ConfigOption
    @Description("""
        Use command line options removed since Docker 1.10.0 (default: `false`).
    """)
    public boolean legacy;

    @ConfigOption
    @Description("""
        Add the specified flags to the volume mounts e.g. `'ro,Z'`.
    """)
    public String mountFlags;

    @ConfigOption
    @Description("""
        The registry from where Docker images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.
    """)
    public String registry;

    @ConfigOption
    @Description("""
        Clean up the container after the execution (default: `true`). See the [Docker documentation](https://docs.docker.com/engine/reference/run/#clean-up---rm) for details.
    """)
    public boolean remove;

    @ConfigOption
    @Description("""
        This attribute can be used to provide any extra command line options supported by the `docker run` command. See the [Docker documentation](https://docs.docker.com/engine/reference/run/) for details.
    """)
    public String runOptions;

    @ConfigOption
    @Description("""
        Executes Docker run command as `sudo` (default: `false`).
    """)
    public boolean sudo;

    @ConfigOption
    @Description("""
        Mounts a path of your choice as the `/tmp` directory in the container. Use the special value `'auto'` to create a temporary directory each time a container is created.
    """)
    public String temp;

    @ConfigOption
    @Description("""
        Allocates a pseudo-tty (default: `false`).
    """)
    public boolean tty;

}
