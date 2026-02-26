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

import java.nio.file.Path;
import java.util.List;

import nextflow.config.schema.ConfigOption;
import nextflow.config.schema.ConfigScope;
import nextflow.script.dsl.Description;
import nextflow.script.types.Duration;

public class CondaConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        Enable Conda execution (default: `false`).
    """)
    public boolean enabled;

    @ConfigOption
    @Description("""
        The path where Conda environments are stored.
    """)
    public Path cacheDir;

    @ConfigOption
    @Description("""
        The Conda channels that can be used to resolve Conda packages.
    """)
    public List<String> channels;

    @ConfigOption
    @Description("""
        Extra command line options to append to the `conda create` command.
    """)
    public String createOptions;

    @ConfigOption
    @Description("""
        The amount of time to wait for the Conda environment to be created before failing (default: `20 min`).
    """)
    public Duration createTimeout;

    @ConfigOption
    @Description("""
        When `true`, use `mamba` instead of `conda` to create the Conda environments.

        [Read more](https://github.com/mamba-org/mamba)
    """)
    public boolean useMamba;

    @ConfigOption
    @Description("""
        When `true`, use `micromamba` instead of `conda` to create the Conda environments.

        [Read more](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
    """)
    public boolean useMicromamba;

}
