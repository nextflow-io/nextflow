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
import nextflow.script.types.Duration;

public class CharliecloudConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        The directory where remote Charliecloud images are stored. When using a computing cluster it must be a shared folder accessible to all compute nodes.
    """)
    public String cacheDir;

    @ConfigOption
    @Description("""
        Enable Charliecloud execution (default: `false`).
    """)
    public boolean enabled;

    @ConfigOption
    @Description("""
        Comma separated list of environment variable names to be included in the container environment.
    """)
    public String envWhitelist;

    @ConfigOption
    @Description("""
        The amount of time the Charliecloud pull can last, exceeding which the process is terminated (default: `20 min`).
    """)
    public Duration pullTimeout;

    @ConfigOption
    @Description("""
        The registry from where images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.
    """)
    public String registry;

    @ConfigOption
    @Description("""
        This attribute can be used to provide any extra command line options supported by the `ch-run` command.
    """)
    public String runOptions;

    @ConfigOption
    @Description("""
        Mounts a path of your choice as the `/tmp` directory in the container. Use the special value `'auto'` to create a temporary directory each time a container is created.
    """)
    public String temp;

    @ConfigOption
    @Description("""
        Create a temporary squashFS container image in the process work directory instead of a folder.
    """)
    public boolean useSquash;

    @ConfigOption
    @Description("""
        Enable `writeFake` with charliecloud. This allows to run containers from storage in writeable mode using overlayfs.
    """)
    public boolean writeFake;

}
