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

import java.util.List;

import nextflow.config.schema.ConfigOption;
import nextflow.config.schema.ConfigScope;
import nextflow.script.dsl.Description;
import nextflow.script.types.MemoryUnit;

public class GoogleBatchConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        The set of allowed locations for VMs to be provisioned (default: no restriction).

        [Read more](https://cloud.google.com/batch/docs/reference/rest/v1/projects.locations.jobs#locationpolicy)
    """)
    public List<String> allowedLocations;

    @ConfigOption
    @Description("""
        The size of the virtual machine boot disk, e.g `50.GB` (default: none).
    """)
    public MemoryUnit bootDiskSize;

    @ConfigOption
    @Description("""
        The minimum CPU Platform, e.g. `'Intel Skylake'` (default: none).

        [Read more](https://cloud.google.com/compute/docs/instances/specify-min-cpu-platform#specifications)
    """)
    public String cpuPlatform;

    @ConfigOption
    @Description("""
        Max number of execution attempts of a job interrupted by a Compute Engine spot reclaim event (default: `5`).
    """)
    public int maxSpotAttempts;

    @ConfigOption
    @Description("""
        The URL of an existing network resource to which the VM will be attached.
    """)
    public String network;

    @ConfigOption
    @Description("""
        The Google service account email to use for the pipeline execution. If not specified, the default Compute Engine service account for the project will be used.

        [Read more](https://www.nextflow.io/docs/latest/google.html#credentials)
    """)
    public String serviceAccountEmail;

    @ConfigOption
    @Description("""
        When `true`, enables the usage of *spot* virtual machines (default: `false`).
    """)
    public boolean spot;

    @ConfigOption
    @Description("""
        The URL of an existing subnetwork resource in the network to which the VM will be attached.
    """)
    public String subnetwork;

    @ConfigOption
    @Description("""
        When `true`, the VM will *not* be provided with a public IP address, and only contain an internal IP.
    """)
    public boolean usePrivateAddress;

}
