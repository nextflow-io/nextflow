/*
 * Copyright 2022, Google Inc.
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

package nextflow.cloud.google.batch.client

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.script.dsl.Description
import nextflow.util.MemoryUnit
/**
 * Model Google Batch config settings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class BatchConfig implements ConfigScope {

    static final private int DEFAULT_MAX_SPOT_ATTEMPTS = 0
    
    static final private List<Integer> DEFAULT_RETRY_LIST = List.of(50001)

    static final private List<String> DEFAULT_GCSFUSE_OPTS = List.<String>of('-o rw', '-implicit-dirs')

    @ConfigOption
    @Description("""
        The set of [allowed locations](https://cloud.google.com/batch/docs/reference/rest/v1/projects.locations.jobs#locationpolicy) for VMs to be provisioned (default: no restriction).
    """)
    final List<String> allowedLocations

    @ConfigOption
    @Description("""
        The list of exit codes that should be automatically retried by Google Batch when `google.batch.maxSpotAttempts` is greater than 0 (default: `[50001]`).

        [Read more](https://cloud.google.com/batch/docs/troubleshooting#reserved-exit-codes)
    """)
    final List<Integer> autoRetryExitCodes

    @ConfigOption
    @Description("""
        The image URI of the virtual machine boot disk, e.g `batch-debian` (default: none).

        [Read more](https://cloud.google.com/batch/docs/vm-os-environment-overview#vm-os-image-options)
    """)
    final String bootDiskImage

    @ConfigOption
    @Description("""
        The size of the virtual machine boot disk, e.g `50.GB` (default: none).
    """)
    final MemoryUnit bootDiskSize

    @ConfigOption
    @Description("""
        The [minimum CPU Platform](https://cloud.google.com/compute/docs/instances/specify-min-cpu-platform#specifications), e.g. `'Intel Skylake'` (default: none).
    """)
    final String cpuPlatform

    @ConfigOption
    @Description("""
        List of custom mount options for `gcsfuse` (default: `['-o rw', '-implicit-dirs']`).
    """)
    final List<String> gcsfuseOptions

    @ConfigOption
    @Description("""
    """)
    final boolean installGpuDrivers

    @ConfigOption
    @Description("""
        Enable the installation of the Ops Agent on Google Batch instances for enhanced monitoring and logging (default: `false`).
    """)
    final boolean installOpsAgent

    @ConfigOption
    @Description("""
        The Google Cloud Storage path where job logs should be stored, e.g. `gs://my-logs-bucket/logs`.
    """)
    final String logsPath

    @ConfigOption
    @Description("""
        Max number of execution attempts of a job interrupted by a Compute Engine Spot reclaim event (default: `0`).
    """)
    final int maxSpotAttempts

    @ConfigOption
    @Description("""
        The URL of an existing network resource to which the VM will be attached.
    """)
    final String network

    @ConfigOption
    @Description("""
        The [network tags](https://cloud.google.com/vpc/docs/add-remove-network-tags) to be applied to the instances created by Google Batch jobs (e.g., `['allow-ssh', 'allow-http']`).
    """)
    final List<String> networkTags

    @ConfigOption
    @Description("""
    """)
    final boolean preemptible

    final BatchRetryConfig retry

    @ConfigOption
    @Description("""
        The Google service account email to use for the pipeline execution. If not specified, the default Compute Engine service account for the project will be used.
    """)
    final String serviceAccountEmail

    @ConfigOption
    @Description("""
        Enable the use of spot virtual machines (default: `false`).
    """)
    final boolean spot

    @ConfigOption
    @Description("""
        The URL of an existing subnetwork resource in the network to which the VM will be attached.
    """)
    final String subnetwork

    @ConfigOption
    @Description("""
        Do not provision public IP addresses for VMs, such that they only have an internal IP address (default: `false`).
    """)
    final boolean usePrivateAddress

    BatchConfig(Map opts) {
        allowedLocations = opts.allowedLocations as List<String> ?: Collections.emptyList()
        autoRetryExitCodes = opts.autoRetryExitCodes as List<Integer> ?: DEFAULT_RETRY_LIST
        bootDiskImage = opts.bootDiskImage
        bootDiskSize = opts.bootDiskSize as MemoryUnit
        cpuPlatform = opts.cpuPlatform
        gcsfuseOptions = opts.gcsfuseOptions as List<String> ?: DEFAULT_GCSFUSE_OPTS
        installGpuDrivers = opts.installGpuDrivers as boolean
        installOpsAgent = opts.installOpsAgent as boolean
        logsPath = opts.logsPath
        maxSpotAttempts = opts.maxSpotAttempts != null ? opts.maxSpotAttempts as int : DEFAULT_MAX_SPOT_ATTEMPTS
        network = opts.network
        networkTags = opts.networkTags as List<String> ?: Collections.emptyList()
        preemptible = opts.preemptible as boolean
        retry = new BatchRetryConfig( opts.retryPolicy as Map ?: Collections.emptyMap() )
        serviceAccountEmail = opts.serviceAccountEmail
        spot = opts.spot as boolean
        subnetwork = opts.subnetwork
        usePrivateAddress = opts.usePrivateAddress as boolean
    }

    BatchRetryConfig getRetryConfig() { retry }

    Path logsPath() {
        return logsPath as Path
    }

}
