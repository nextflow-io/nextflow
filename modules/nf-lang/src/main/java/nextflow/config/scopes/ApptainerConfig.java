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

import groovy.transform.CompileStatic;
import nextflow.config.schema.ConfigOption;
import nextflow.config.schema.ConfigScope;
import nextflow.script.dsl.Description;
import nextflow.script.types.Duration;

public class ApptainerConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        When `true` Nextflow automatically mounts host paths in the executed container. It requires the `user bind control` feature to be enabled in your Apptainer installation (default: `true`).
    """)
    public boolean autoMounts;

    @ConfigOption
    @Description("""
        The directory where remote Apptainer images are stored. When using a computing cluster it must be a shared folder accessible to all compute nodes.
    """)
    public String cacheDir;

    @ConfigOption
    @Description("""
        Enable Apptainer execution (default: `false`).
    """)
    public boolean enabled;

    @ConfigOption
    @Description("""
        This attribute can be used to provide any option supported by the Apptainer engine i.e. `apptainer [OPTIONS]`.
    """)
    public String engineOptions;

    @ConfigOption
    @Description("""
        Comma separated list of environment variable names to be included in the container environment.
    """)
    public String envWhitelist;

    @ConfigOption
    @Description("""
        Directory where remote Apptainer images are retrieved. When using a computing cluster it must be a shared folder accessible to all compute nodes.
    """)
    public String libraryDir;

    @ConfigOption
    @Description("""
        Pull the Apptainer image with http protocol (default: `false`).
    """)
    public boolean noHttps;

    @ConfigOption
    @Description("""
        When enabled, OCI (and Docker) container images are pulled and converted to the SIF format by the Apptainer run command, instead of Nextflow (default: `false`).
    """)
    public boolean ociAutoPull;

    @ConfigOption
    @Description("""
        The amount of time the Apptainer pull can last, exceeding which the process is terminated (default: `20 min`).
    """)
    public Duration pullTimeout;

    @ConfigOption
    @Description("""
        The registry from where Docker images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.
    """)
    public String registry;

    @ConfigOption
    @Description("""
        This attribute can be used to provide any extra command line options supported by `apptainer exec`.
    """)
    public String runOptions;

}
