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

import nextflow.config.schema.ConfigOption;
import nextflow.config.schema.ConfigScope;
import nextflow.script.dsl.Description;
import nextflow.script.types.Duration;

public class SpackConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        The path where Spack environments are stored.
    """)
    public Path cacheDir;

    @ConfigOption
    @Description("""
        Enables checksum verification for source tarballs (default: `true`).
    """)
    public boolean checksum;

    @ConfigOption
    @Description("""
        The amount of time to wait for the Spack environment to be created before failing (default: `60 min`).
    """)
    public Duration createTimeout;

    @ConfigOption
    @Description("""
        The maximum number of parallel package builds (default: the number of available CPUs).
    """)
    public int parallelBuilds;

}
