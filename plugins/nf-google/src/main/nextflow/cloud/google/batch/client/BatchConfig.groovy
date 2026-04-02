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

package nextflow.cloud.google.batch.client

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.CloudTransferOptions
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.script.dsl.Description
import nextflow.util.Duration
import nextflow.util.MemoryUnit
/**
 * Model Google Batch config settings.
 * <p>
 * <strong>Copy transports</strong> ({@link #stageInCopyTransport}, {@link #stageOutCopyTransport}): optional values
 * {@code posix}, {@code gcloud}, {@code gsutil}. When unset or only {@code posix}, staging uses
 * {@link nextflow.executor.SimpleFileCopyStrategy} with gcsfuse paths under {@code /mnt/disks}. When {@code gcloud}
 * or {@code gsutil} is set for a direction, {@link nextflow.cloud.google.batch.GoogleBatchFileCopyStrategy} loads and
 * generated bash uses {@code gcloud storage} / {@code gsutil} with fallbacks to POSIX copy from the mount; order is
 * chosen per transport (see {@link nextflow.cloud.google.batch.GoogleBatchBashLib}). CLI stage-in applies when
 * {@code process stageInMode} is {@code copy}; CLI stage-out when effective {@code stageOutMode} is {@code copy}.
 * {@code move} / {@code rsync} / etc. keep existing POSIX behaviour. Parallelism and retries use
 * {@link #maxParallelTransfers}, {@link #maxTransferAttempts}, {@link #delayBetweenAttempts} with {@code nxf_parallel}
 * / {@code nxf_cp_retry} in the task script.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class BatchConfig implements ConfigScope {

    /**
     * Copy via POSIX {@code cp}/links from gcsfuse paths under {@code /mnt/disks/...}.
     */
    static final String COPY_TRANSPORT_POSIX = 'posix'

    /**
     * Copy using {@code gcloud storage} first, then {@code gsutil}, then POSIX from the mount (see {@link nextflow.cloud.google.batch.GoogleBatchBashLib}).
     */
    static final String COPY_TRANSPORT_GCLOUD = 'gcloud'

    /**
     * Copy using {@code gsutil} first, then {@code gcloud storage}, then POSIX from the mount.
     */
    static final String COPY_TRANSPORT_GSUTIL = 'gsutil'

    static final private List<String> VALID_COPY_TRANSPORTS = List.of(COPY_TRANSPORT_POSIX, COPY_TRANSPORT_GCLOUD, COPY_TRANSPORT_GSUTIL)

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
        `gcloud` executable path (default: `gcloud` on `PATH`).
    """)
    final String gcloudCli

    @ConfigOption
    @Description("""
        `gsutil` executable path (default: `gsutil` on `PATH`).
    """)
    final String gsutilCli

    @ConfigOption
    @Description("""
        Max parallel CLI object copies (default: same as other cloud executors).
    """)
    final int maxParallelTransfers

    @ConfigOption
    @Description("""
        Max retries per CLI copy (default: same as other cloud executors).
    """)
    final int maxTransferAttempts

    @ConfigOption
    @Description("""
        Delay between CLI copy retries (default: same as other cloud executors).
    """)
    final Duration delayBetweenAttempts

    @ConfigOption
    @Description("""
        Stage-in copy transport: `posix`, `gcloud`, or `gsutil`.
    """)
    final String stageInCopyTransport

    @ConfigOption
    @Description("""
        Stage-out copy transport: `posix`, `gcloud`, or `gsutil`.
    """)
    final String stageOutCopyTransport

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
        gcloudCli = opts.gcloudCli as String
        gsutilCli = opts.gsutilCli as String
        maxParallelTransfers = opts.maxParallelTransfers != null ? opts.maxParallelTransfers as int : CloudTransferOptions.MAX_TRANSFER
        maxTransferAttempts = opts.maxTransferAttempts != null ? opts.maxTransferAttempts as int : CloudTransferOptions.MAX_TRANSFER_ATTEMPTS
        delayBetweenAttempts = opts.delayBetweenAttempts ? opts.delayBetweenAttempts as Duration : CloudTransferOptions.DEFAULT_DELAY_BETWEEN_ATTEMPTS
        stageInCopyTransport = opts.stageInCopyTransport as String
        stageOutCopyTransport = opts.stageOutCopyTransport as String
        for (String t in [stageInCopyTransport, stageOutCopyTransport]) {
            if (t && t !in VALID_COPY_TRANSPORTS)
                throw new IllegalArgumentException("Invalid google.batch copy transport: '$t' — valid values are: ${VALID_COPY_TRANSPORTS.join(', ')}")
        }
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

    /**
     * Load {@link nextflow.cloud.google.batch.GoogleBatchFileCopyStrategy} when any transport requests CLI object-storage copy.
     */
    boolean usesGoogleBatchStaging() {
        return [stageInCopyTransport, stageOutCopyTransport].any { it in [COPY_TRANSPORT_GCLOUD, COPY_TRANSPORT_GSUTIL] }
    }

    BatchRetryConfig getRetryConfig() { retry }

    Path logsPath() {
        return logsPath as Path
    }

}
