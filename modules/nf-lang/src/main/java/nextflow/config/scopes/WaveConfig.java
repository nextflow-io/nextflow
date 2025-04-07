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

public class WaveConfig implements ConfigScope {

    public WaveBuildConfig build;

    @ConfigOption
    @Description("""
        Enable the use of Wave containers.
    """)
    public boolean enabled;

    @ConfigOption
    @Description("""
        The Wave service endpoint (default: `https://wave.seqera.io`).
    """)
    public String endpoint;

    @ConfigOption
    @Description("""
        Enables Wave container freezing. Wave will provision a non-ephemeral container image that will be pushed to a container repository of your choice.

        See also: `wave.build.repository` and `wave.build.cacheRepository`
    """)
    public boolean freeze;

    public WaveHttpConfig http;

    @ConfigOption
    @Description("""
        Enables Wave container mirroring.
    """)
    public boolean mirror;

    public WaveRetryConfig retryPolicy;

    public WaveScanConfig scan;

    @ConfigOption
    @Description("""
        The strategy to be used when resolving multiple Wave container requirements (default: `'container,dockerfile,conda,spack'`).
    """)
    public String strategy;

}
