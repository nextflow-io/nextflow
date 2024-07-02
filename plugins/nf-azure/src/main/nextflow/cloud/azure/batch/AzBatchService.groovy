/*
 * Copyright 2021, Microsoft Corp
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

package nextflow.cloud.azure.batch

import java.math.RoundingMode
import java.nio.file.Path
import java.time.Duration
import java.time.Instant
import java.time.temporal.ChronoUnit
import java.util.concurrent.TimeoutException
import java.util.function.Predicate

import com.azure.compute.batch.BatchClient
import com.azure.compute.batch.BatchClientBuilder
import com.azure.compute.batch.models.AutoUserScope
import com.azure.compute.batch.models.AutoUserSpecification
import com.azure.compute.batch.models.AzureFileShareConfiguration
import com.azure.compute.batch.models.BatchJobCreateContent
import com.azure.compute.batch.models.BatchJobUpdateContent
import com.azure.compute.batch.models.BatchNodeFillType
import com.azure.compute.batch.models.BatchPool
import com.azure.compute.batch.models.BatchPoolCreateContent
import com.azure.compute.batch.models.BatchPoolInfo
import com.azure.compute.batch.models.BatchPoolState
import com.azure.compute.batch.models.BatchStartTask
import com.azure.compute.batch.models.BatchSupportedImage
import com.azure.compute.batch.models.BatchTask
import com.azure.compute.batch.models.BatchTaskConstraints
import com.azure.compute.batch.models.BatchTaskContainerSettings
import com.azure.compute.batch.models.BatchTaskCreateContent
import com.azure.compute.batch.models.BatchTaskSchedulingPolicy
import com.azure.compute.batch.models.ContainerConfiguration
import com.azure.compute.batch.models.ContainerRegistryReference
import com.azure.compute.batch.models.ContainerType
import com.azure.compute.batch.models.ElevationLevel
import com.azure.compute.batch.models.MetadataItem
import com.azure.compute.batch.models.MountConfiguration
import com.azure.compute.batch.models.NetworkConfiguration
import com.azure.compute.batch.models.OnAllBatchTasksComplete
import com.azure.compute.batch.models.OutputFile
import com.azure.compute.batch.models.OutputFileBlobContainerDestination
import com.azure.compute.batch.models.OutputFileDestination
import com.azure.compute.batch.models.OutputFileUploadCondition
import com.azure.compute.batch.models.OutputFileUploadConfig
import com.azure.compute.batch.models.ResourceFile
import com.azure.compute.batch.models.UserIdentity
import com.azure.compute.batch.models.VirtualMachineConfiguration
import com.azure.core.credential.AzureNamedKeyCredential
import com.azure.core.credential.TokenCredential
import com.azure.core.exception.ClientAuthenticationException
import com.azure.core.exception.HttpResponseException
import com.azure.core.exception.ResourceNotFoundException
import com.azure.core.http.rest.PagedIterable
import com.azure.identity.ClientSecretCredentialBuilder
import com.azure.identity.ManagedIdentityCredentialBuilder
import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedSupplier
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
import nextflow.cloud.azure.config.AzConfig
import nextflow.cloud.azure.config.AzFileShareOpts
import nextflow.cloud.azure.config.AzPoolOpts
import nextflow.cloud.azure.config.AzStartTaskOpts
import nextflow.cloud.azure.config.CopyToolInstallMode
import nextflow.cloud.azure.nio.AzPath
import nextflow.cloud.types.CloudMachineInfo
import nextflow.cloud.types.PriceModel
import nextflow.fusion.FusionHelper
import nextflow.fusion.FusionScriptLauncher
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.util.CacheHelper
import nextflow.util.MemoryUnit
import nextflow.util.MustacheTemplateEngine
import nextflow.util.Rnd
/**
 * Implements Azure Batch operations for Nextflow executor
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AzBatchService implements Closeable {

    static private final String AZCOPY_URL = 'https://github.com/nextflow-io/azcopy-tool/raw/linux_amd64_10.8.0/azcopy'

    static private final long _1GB = 1 << 30

    static final private Map<String,AzVmPoolSpec> allPools = new HashMap<>(50)

    AzConfig config

    Map<TaskProcessor,String> allJobIds = new HashMap<>(50)

    AzBatchService(AzBatchExecutor executor) {
        assert executor
        this.config = executor.config
    }

    protected AzVmPoolSpec getPoolSpec(String poolId) {
        return allPools.get(poolId)
    }

    @Memoized
    protected BatchClient getClient() {
        createBatchClient()
    }

    /**
     * @return The available location codes
     */
    @Memoized
    List<String> listLocationNames() {
        // this has been generated with `az account list-locations` cli tool
        // it needs to be replaced by a proper API call
        final json = AzBatchService.class.getResourceAsStream('/nextflow/cloud/azure/az-locations.json')
        if( !json )
            throw new IllegalArgumentException("Unable to fetch Azure locations")
        final locs = (new JsonSlurper().parse(json) as List<Map>)
        return locs.collect {  (String)it.name }
    }

    @Memoized
    private List<Map> listAllVms(String location) {
        if( !location )
            throw new IllegalArgumentException("Missing Azure location parameter")
        final json = AzBatchService.class.getResourceAsStream("/nextflow/cloud/azure/vm-list-size-${location}.json")
        if( !json ) {
            log.warn "Unable to find Azure VM names for location: $location"
            return Collections.emptyList()
        }

        return (List<Map>) new JsonSlurper().parse(json)
    }

    @Memoized
    List<String> listVmNames(String location) {
        return listAllVms(location).collect { it.name as String }
    }

    AzVmType guessBestVm(String location, int cpus, MemoryUnit mem, String family) {
        log.debug "[AZURE BATCH] guessing best VM given location=$location; cpus=$cpus; mem=$mem; family=$family"
        if( !family.contains('*') && !family.contains('?') )
            return findBestVm(location, cpus, mem, family)

        // well this is a quite heuristic tentative to find a bigger instance to accommodate more tasks
        AzVmType result=null
        if( cpus<=4 ) {
            result = findBestVm(location, cpus*4, mem!=null ? mem*4 : null, family)
            if( !result )
                result = findBestVm(location, cpus*2, mem!=null ? mem*2 : null, family)
        }
        else if( cpus <=8 ) {
            result = findBestVm(location, cpus*2, mem!=null ? mem*2 : null, family)
        }
        if( !result )
            result = findBestVm(location, cpus, mem, family)
        return result
    }

    /**
     * Find the best VM matching the specified criteria
     *
     * @param location The azure location
     * @param cpus The requested number of cpus
     * @param mem The requested amount of memory. Can be null if no mem has been specified
     * @param allFamilies Comma separate list of Azure VM machine types, each value can also contain wildcard characters ie. `*` and `?`
     * @return The `AzVmType` instance that best accommodate the resource requirement
     */
    AzVmType findBestVm(String location, int cpus, MemoryUnit mem, String allFamilies) {
        def all = listAllVms(location)
        def scores = new TreeMap<Double,String>()
        def list = allFamilies ? allFamilies.tokenize(',') : ['']
        for( String family : list ) {
            for( Map entry : all ) {
                if( !matchType(family, entry.name as String) )
                    continue
                def score = computeScore(cpus, mem, entry)
                if( score != null )
                    scores.put(score, entry.name as String)
            }
        }

        return scores ? getVmType(location, scores.firstEntry().value) : null
    }

    protected boolean matchType(String family, String vmType) {
        if( !family )
            return true
        if( family.contains('*') )
            family = family.toLowerCase().replaceAll(/\*/,'.*')
        else if( family.contains('?') )
            family = family.toLowerCase().replaceAll(/\?/,'.{1}')

        return vmType =~ /(?i)^${family}$/
    }

    protected Double computeScore(int cpus, MemoryUnit mem, Map entry) {
        def vmCores = entry.numberOfCores as int
        double vmMemGb = (entry.memoryInMB as int) /1024

        if( cpus > vmCores ) {
            return null
        }

        int cpusDelta = cpus-vmCores
        double score = cpusDelta * cpusDelta
        if( mem && vmMemGb ) {
            double memGb = mem.toMega()/1024
            if( memGb > vmMemGb )
                return null
            double memDelta = memGb - vmMemGb
            score += memDelta*memDelta
        }

        return Math.sqrt(score)
    }

    @Memoized
    AzVmType getVmType(String location, String vmName) {
        def vm = listAllVms(location).find { vmName.equalsIgnoreCase(it.name?.toString()) }
        if( !vm )
            throw new IllegalArgumentException("Unable to find size for VM name '$vmName' and location '$location'")

        new AzVmType(vm)
    }

    protected int computeSlots(int cpus, MemoryUnit mem, int vmCpus, MemoryUnit vmMem) {
        //  cpus requested should not exceed max cpus avail
        final cpuSlots = Math.min(cpus, vmCpus) as int
        if( !mem || !vmMem )
            return cpuSlots
        //  mem requested should not exceed max mem avail
        final vmMemGb = vmMem.mega /_1GB as float
        final memGb = mem.mega /_1GB as float
        final mem0 = Math.min(memGb, vmMemGb)
        return Math.max(cpuSlots, memSlots(mem0, vmMemGb, vmCpus))
    }

    protected int computeSlots(TaskRun task, AzVmPoolSpec pool) {
        computeSlots(
                task.config.getCpus(),
                task.config.getMemory(),
                pool.vmType.numberOfCores,
                pool.vmType.memory )
    }


    protected int memSlots(float memGb, float vmMemGb, int vmCpus) {
        BigDecimal result = memGb / (vmMemGb / vmCpus)
        result.setScale(0, RoundingMode.UP).intValue()
    }


    protected AzureNamedKeyCredential createBatchCredentialsWithKey() {
        log.debug "[AZURE BATCH] Creating Azure Batch client using shared key creddentials"

        if( !config.batch().endpoint )
            throw new IllegalArgumentException("Missing Azure Batch endpoint -- Specify it in the nextflow.config file using the setting 'azure.batch.endpoint'")
        if( !config.batch().accountName )
            throw new IllegalArgumentException("Missing Azure Batch account name -- Specify it in the nextflow.config file using the setting 'azure.batch.accountName'")
        if( !config.batch().accountKey )
            throw new IllegalArgumentException("Missing Azure Batch account key -- Specify it in the nextflow.config file using the setting 'azure.batch.accountKey'")

        return new AzureNamedKeyCredential(config.batch().accountName, config.batch().accountKey)
    }

    protected TokenCredential createBatchCredentialsWithServicePrincipal() {
        log.debug "[AZURE BATCH] Creating Azure Batch client using Service Principal credentials"

        return new ClientSecretCredentialBuilder()
                .tenantId(config.activeDirectory().tenantId)
                .clientId(config.activeDirectory().servicePrincipalId)
                .clientSecret(config.activeDirectory().servicePrincipalSecret)
                .build()
    }

    protected TokenCredential createBatchCredentialsWithManagedIdentity() {
        final clientId = config.managedIdentity().clientId
        final credential = new ManagedIdentityCredentialBuilder()
        if (clientId) {
            log.debug "[AZURE BATCH] Creating Azure Batch client using Managed Identity credentials - clientId: ${clientId}"
            credential.clientId(clientId)
        }
        else {
            log.debug '[AZURE BATCH] Creating Azure Batch client using Managed Identity credentials - not clientId provided'
        }
        return credential.build()
    }

    protected BatchClient createBatchClient() {
        log.debug "[AZURE BATCH] Executor options=${config.batch()}"

        final builder = new BatchClientBuilder()
        if( config.managedIdentity().isConfigured() )
            builder.credential( createBatchCredentialsWithManagedIdentity() )
        else if( config.activeDirectory().isConfigured() )
            builder.credential( createBatchCredentialsWithServicePrincipal() )
        else if( config.batch().endpoint || config.batch().accountKey || config.batch().accountName )
            builder.credential( createBatchCredentialsWithKey() )

        if( config.batch().endpoint )
            builder.endpoint(config.batch().endpoint)

        return builder.buildClient()
    }

    AzTaskKey submitTask(TaskRun task) {
        final poolId = getOrCreatePool(task)
        final jobId = getOrCreateJob(poolId, task)
        runTask(poolId, jobId, task)
    }

    BatchTask getTask(AzTaskKey key) {
        apply(() -> client.getTask(key.jobId, key.taskId))
    }

    void terminate(AzTaskKey key) {
        apply(() -> client.terminateTask(key.jobId, key.taskId))
    }

    CloudMachineInfo machineInfo(AzTaskKey key) {
        if( !key || !key.jobId )
            throw new IllegalArgumentException("Missing Azure Batch job id")
        final job = apply(() -> client.getJob(key.jobId))
        final poolId = job.poolInfo.poolId
        final AzVmPoolSpec spec = allPools.get(poolId)
        if( !spec )
            throw new IllegalStateException("Missing VM pool spec for pool id: $poolId")
        final type = spec.vmType.name
        if( !type )
            throw new IllegalStateException("Missing VM type for pool id: $poolId")
        new CloudMachineInfo(type: type, priceModel: PriceModel.standard, zone: config.batch().location)
    }

    synchronized String getOrCreateJob(String poolId, TaskRun task) {
        final mapKey = task.processor
        if( allJobIds.containsKey(mapKey)) {
            return allJobIds[mapKey]
        }
        // create a batch job
        final jobId = makeJobId(task)
        final content = new BatchJobCreateContent(jobId, new BatchPoolInfo(poolId: poolId))
        apply(() -> client.createJob(content))
        // add to the map
        allJobIds[mapKey] = jobId
        return jobId
    }

    String makeJobId(TaskRun task) {
        final name = task
                .processor
                .name
                .trim()
                .replaceAll(/[^a-zA-Z0-9-_]+/,'_')

        final String key = "job-${Rnd.hex()}-${name}"
        // Azure batch job max len is 64 characters, however we keep it a bit shorter
        // because the jobId + taskId composition must be less then 100
        final MAX_LEN = 62i
        return key.size()>MAX_LEN ? key.substring(0,MAX_LEN) : key
    }

    protected BatchTaskCreateContent createTask(String poolId, String jobId, TaskRun task) {
        assert poolId, 'Missing Azure Batch poolId argument'
        assert jobId, 'Missing Azure Batch jobId argument'
        assert task, 'Missing Azure Batch task argument'

        final sas = config.storage().sasToken
        if( !sas )
            throw new IllegalArgumentException("Missing Azure Blob storage SAS token")

        final container = task.getContainer()
        if( !container )
            throw new IllegalArgumentException("Missing container image for process: $task.name")
        final taskId = "nf-${task.hash.toString()}"
        // get the pool config
        final pool = getPoolSpec(poolId)
        if( !pool )
            throw new IllegalStateException("Missing Azure Batch pool spec with id: $poolId")
        // container settings
        // mount host certificates otherwise `azcopy` fails
        def opts = "-v /etc/ssl/certs:/etc/ssl/certs:ro -v /etc/pki:/etc/pki:ro "
        // shared volume mounts
        final shares = getShareVolumeMounts(pool)
        if( shares )
            opts += "${shares.join(' ')} "
        // custom container settings
        if( task.config.getContainerOptions() )
            opts += "${task.config.getContainerOptions()} "
        // fusion environment settings
        final fusionEnabled = FusionHelper.isFusionEnabled((Session)Global.session)
        final launcher = fusionEnabled ? FusionScriptLauncher.create(task.toTaskBean(), 'az') : null
        if( fusionEnabled ) {
            opts += "--privileged "
            for( Map.Entry<String,String> it : launcher.fusionEnv() ) {
                opts += "-e $it.key=$it.value "
            }
        }
        // config overall container settings
        final containerOpts = new BatchTaskContainerSettings(container)
                .setContainerRunOptions(opts)
        // submit command line
        final String cmd = fusionEnabled
                ? launcher.fusionSubmitCli(task).join(' ')
                : "sh -c 'bash ${TaskRun.CMD_RUN} 2>&1 | tee ${TaskRun.CMD_LOG}'"
        // cpus and memory
        final slots = computeSlots(task, pool)
        // max wall time
        final constraints = new BatchTaskConstraints()
        if( task.config.getTime() )
            constraints.setMaxWallClockTime( Duration.of(task.config.getTime().toMillis(), ChronoUnit.MILLIS) )

        log.trace "[AZURE BATCH] Submitting task: $taskId, cpus=${task.config.getCpus()}, mem=${task.config.getMemory()?:'-'}, slots: $slots"

        return new BatchTaskCreateContent(taskId, cmd)
                .setUserIdentity(userIdentity(pool.opts.privileged, pool.opts.runAs, AutoUserScope.TASK))
                .setContainerSettings(containerOpts)
                .setResourceFiles(resourceFileUrls(task, sas))
                .setOutputFiles(outputFileUrls(task, sas))
                .setRequiredSlots(slots)
                .setConstraints(constraints)
                

    }

    AzTaskKey runTask(String poolId, String jobId, TaskRun task) {
        final taskToAdd = createTask(poolId, jobId, task)
        apply(() -> client.createTask(jobId, taskToAdd))
        return new AzTaskKey(jobId, taskToAdd.getId())
    }

    protected List<String> getShareVolumeMounts(AzVmPoolSpec spec) {
        assert spec!=null
        final shares = config.storage().getFileShares()
        if( !shares )
            return Collections.<String>emptyList()

        // get the file share root
        final mountPath = spec.opts?.getFileShareRootPath()
        if ( !mountPath )
            throw new IllegalArgumentException("Missing FileShareRootPath for pool: ${spec.poolId}")

        final result = new ArrayList(shares.size())
        for( Map.Entry<String, AzFileShareOpts> it : shares ) {
            result.add("-v ${mountPath}/${it.key}:${it.value.mountPath}:rw")
        }
        return result
    }

    protected List<ResourceFile> resourceFileUrls(TaskRun task, String sas) {
        final cmdRun = (AzPath) task.workDir.resolve(TaskRun.CMD_RUN)
        final cmdScript = (AzPath) task.workDir.resolve(TaskRun.CMD_SCRIPT)

        final resFiles = new ArrayList(10)

        if( config.batch().copyToolInstallMode == CopyToolInstallMode.task ) {
            log.trace "[AZURE BATCH] installing azcopy as task resource"
            resFiles << new ResourceFile()
                    .setHttpUrl(AZCOPY_URL)
                    .setFilePath('.nextflow-bin/azcopy')
        }

        resFiles << new ResourceFile()
                .setHttpUrl(AzHelper.toHttpUrl(cmdRun, sas))
                .setFilePath(TaskRun.CMD_RUN)

        resFiles << new ResourceFile()
                .setHttpUrl(AzHelper.toHttpUrl(cmdScript, sas))
                .setFilePath(TaskRun.CMD_SCRIPT)

        if( task.stdin ) {
            resFiles << new ResourceFile()
                    .setHttpUrl(AzHelper.toHttpUrl(cmdScript, sas))
                    .setFilePath(TaskRun.CMD_INFILE)
        }

        return resFiles
    }

    protected List<OutputFile> outputFileUrls(TaskRun task, String sas) {
        List<OutputFile> result = new ArrayList<>(20)
        result << destFile(TaskRun.CMD_EXIT, task.workDir, sas)
        result << destFile(TaskRun.CMD_LOG, task.workDir, sas)
        return result
    }

    protected OutputFile destFile(String localPath, Path targetDir, String sas) {
        log.debug "Task output path: $localPath -> ${targetDir.toUriString()}"
        def target = targetDir.resolve(localPath)
        final dest = new OutputFileBlobContainerDestination(AzHelper.toContainerUrl(targetDir,sas))
                .setPath(target.subpath(1,target.nameCount).toString())

        return new OutputFile(localPath, new OutputFileDestination().setContainer(dest), new OutputFileUploadConfig(OutputFileUploadCondition.TASK_COMPLETION))
    }

    protected BatchSupportedImage getImage(AzPoolOpts opts) {
        PagedIterable<BatchSupportedImage> images = apply(() -> client.listSupportedImages())

        for (BatchSupportedImage it : images) {
            if( !it.nodeAgentSkuId.equalsIgnoreCase(opts.sku) )
                continue
            if( it.osType != opts.osType )
                continue
            if( it.verificationType != opts.verification )
                continue
            if( !it.imageReference.publisher.equalsIgnoreCase(opts.publisher) )
                continue
            if( it.imageReference.offer.equalsIgnoreCase(opts.offer) )
                return it
        }

        throw new IllegalStateException("Cannot find a matching VM image with publisher=$opts.publisher; offer=$opts.offer; OS type=$opts.osType; verification type=$opts.verification")
    }

    protected AzVmPoolSpec specFromPoolConfig(String poolId) {

        def opts = config.batch().pool(poolId)
        def name = opts ? opts.vmType : getPool(poolId)?.vmSize
        if( !name )
            throw new IllegalArgumentException("Cannot find Azure Batch config for pool: $poolId")
        if( !opts )
            opts = new AzPoolOpts(vmType: name)

        def type = getVmType(config.batch().location, name)
        if( !type )
            throw new IllegalArgumentException("Cannot find Azure Batch VM type '$poolId' - Check pool definition $poolId in the Nextflow config file")

        return new AzVmPoolSpec(poolId: poolId, vmType: type, opts: opts)
    }

    protected AzVmPoolSpec specFromAutoPool(TaskRun task) {
        // define a stable pool name
        final loc = config.batch().location
        if( !loc )
            throw new IllegalArgumentException("Missing Azure Batch location")

        final opts = config.batch().autoPoolOpts()
        final mem = task.config.getMemory()
        final cpus = task.config.getCpus()
        final type = task.config.getMachineType() ?: opts.vmType
        if( !type )
            throw new IllegalArgumentException("Missing Azure Batch VM type for task '${task.name}'")

        final vmType = guessBestVm(loc, cpus, mem, type)
        if( !vmType ) {
            def msg = "Cannot find a VM for task '${task.name}' matching these requirements: type=$type, cpus=${cpus}, mem=${mem?:'-'}, location=${loc}"
            throw new IllegalArgumentException(msg)
        }

        final metadata = task.config.getResourceLabels()
        final key = CacheHelper.hasher([vmType.name, opts, metadata]).hash().toString()
        final poolId = "nf-pool-$key-$vmType.name"
        return new AzVmPoolSpec(poolId: poolId, vmType: vmType, opts: opts, metadata: metadata)
    }

    protected void checkPool(BatchPool pool, AzVmPoolSpec spec) {
        if( pool.state != BatchPoolState.ACTIVE ) {
            throw new IllegalStateException("Azure Batch pool '${pool.id}' not in active state")
        }
        else if ( pool.resizeErrors && pool.currentDedicatedNodes==0 ) {
            throw new IllegalStateException("Azure Batch pool '${pool.id}' has resize errors")
        }
        if( pool.taskSlotsPerNode != spec.vmType.numberOfCores ) {
            throw new IllegalStateException("Azure Batch pool '${pool.id}' slots per node does not match the VM num cores (slots: ${pool.taskSlotsPerNode}, cores: ${spec.vmType.numberOfCores})")
        }
    }

    protected void checkPoolId(String poolId) {
        if( !poolId.matches(/^[\w\-]+$/) )
            throw new IllegalArgumentException("Invalid Azure Batch pool Id '$poolId' - It can only contain alphanumeric, hyphen and underscore characters")
    }

    protected AzVmPoolSpec specForTask(TaskRun task) {
        String poolId = null
        if( !config.batch().autoPoolMode ) {
            // the process queue is used as poolId
            poolId = task.config.queue as String
            if( !poolId ) {
                throw new IllegalArgumentException("No Azure Batch pool was specified for task '${task.name}' - Either specify the pool name using the 'queue' directive or enable the 'autoPoolMode' option")
            }
            // sanity check
            checkPoolId(poolId)
            // check if cached
            if( allPools.containsKey(poolId) ) {
                return allPools.get(poolId)
            }
        }

        return poolId
                ? specFromPoolConfig(poolId)
                : specFromAutoPool(task)

    }

    synchronized String getOrCreatePool(TaskRun task) {

        final spec = specForTask(task)
        if( spec && allPools.containsKey(spec.poolId) )
            return spec.poolId

        // check existence and create if needed
        log.debug "[AZURE BATCH] Checking VM pool id=$spec.poolId; size=$spec.vmType"
        def pool = getPool(spec.poolId)
        if( !pool ) {
            if( config.batch().canCreatePool() ) {
                createPool(spec)
            }
            else {
                throw new IllegalArgumentException("Can't find Azure Batch pool '$spec.poolId' - Make sure it exists or set `allowPoolCreation=true` in the nextflow config file")
            }
        }
        else {
            checkPool(pool, spec)
        }

        // add to the list of pool ids
        allPools[spec.poolId] = spec

        return spec.poolId
    }


    protected BatchPool getPool(String poolId) {
        try {
            return apply(() -> client.getPool(poolId))
        }
        catch (IllegalArgumentException | ClientAuthenticationException e) {
            throw e
        }
        catch (ResourceNotFoundException ignored) {
            return null
        }
    }

    protected VirtualMachineConfiguration poolVmConfig(AzPoolOpts opts) {
        /**
         * A container configuration must be provided for a task to run in a specific container.
         * Such container can be pre-fetched on VM creation or when running the task
         *
         * https://github.com/MicrosoftDocs/azure-docs/blob/master/articles/batch/batch-docker-container-workloads.md#:~:text=Run%20container%20applications%20on%20Azure,compatible%20containers%20on%20the%20nodes.
         */
        final containerConfig = new ContainerConfiguration(ContainerType.DOCKER_COMPATIBLE)
        final registryOpts = config.registry()

        if( registryOpts && registryOpts.isConfigured() ) {
            final containerRegistries = new ArrayList<ContainerRegistryReference>(1)
            containerRegistries << new ContainerRegistryReference()
                    .setRegistryServer(registryOpts.server)
                    .setUsername(registryOpts.userName)
                    .setPassword(registryOpts.password)
            containerConfig.setContainerRegistries(containerRegistries)
            log.debug "[AZURE BATCH] Connecting Azure Batch pool to Container Registry '$registryOpts.server'"
        }

        final image = getImage(opts)

        new VirtualMachineConfiguration(image.imageReference, image.nodeAgentSkuId)
                .setContainerConfiguration(containerConfig)
    }


    protected BatchStartTask createStartTask(AzStartTaskOpts opts) {
        log.trace "AzStartTaskOpts: ${opts}"
        final startCmd = new ArrayList<String>(5)
        final resourceFiles = new ArrayList<ResourceFile>()

        // If enabled, append azcopy installer to start task command
        if( config.batch().getCopyToolInstallMode() == CopyToolInstallMode.node ) {
            startCmd << 'bash -c "chmod +x azcopy && mkdir \$AZ_BATCH_NODE_SHARED_DIR/bin/ && cp azcopy \$AZ_BATCH_NODE_SHARED_DIR/bin/"'
            resourceFiles << new ResourceFile()
                .setHttpUrl(AZCOPY_URL)
                .setFilePath('azcopy')
        }

        // Get any custom start task command
        if ( opts.script ) {
            startCmd << "bash -c '${opts.script.replace(/'/,/''/)}'".toString()
        }

        // If there is no start task contents we return a null to indicate no start task
        if( !startCmd ) {
            return null
        }

        // otherwise return a StartTask object with the start task command and resource files
        return new BatchStartTask(startCmd.join('; '))
            .setResourceFiles(resourceFiles)
            .setUserIdentity(userIdentity(opts.privileged, null, AutoUserScope.POOL))
    }

    protected void createPool(AzVmPoolSpec spec) {

        final poolParams = new BatchPoolCreateContent(spec.poolId, spec.vmType.name)
                .setVirtualMachineConfiguration(poolVmConfig(spec.opts))
                // same as the number of cores
                // https://docs.microsoft.com/en-us/azure/batch/batch-parallel-node-tasks
                .setTaskSlotsPerNode(spec.vmType.numberOfCores)

        final startTask = createStartTask(spec.opts.startTask)
        if( startTask ) {
            poolParams .setStartTask(startTask)
        }

        // resource labels
        if( spec.metadata ) {
            final metadata = spec.metadata.collect { name, value ->
                new MetadataItem(name, value)
            }
            poolParams.setMetadata(metadata)
        }

        // virtual network
        if( spec.opts.virtualNetwork )
            poolParams.setNetworkConfiguration( new NetworkConfiguration().setSubnetId(spec.opts.virtualNetwork) )

        // scheduling policy
        if( spec.opts.schedulePolicy ) {
            final pol = BatchNodeFillType.fromString(spec.opts.schedulePolicy)
            if( !pol ) throw new IllegalArgumentException("Unknown Azure Batch scheduling policy: ${spec.opts.schedulePolicy}")
            poolParams.setTaskSchedulingPolicy( new BatchTaskSchedulingPolicy(pol) )
        }

        // mount points
        if ( config.storage().fileShares ) {
            List<MountConfiguration> mountConfigs = new ArrayList(config.storage().fileShares.size())
            config.storage().fileShares.each {
                if (it.key) {
                    final String accountName = config.storage().accountName
                    final endpoint = "https://${config.storage().accountName}.file.core.windows.net/${it.key}" as String
                    final accountKey = config.storage().accountKey
                    final shareConfig = new AzureFileShareConfiguration( accountName, endpoint, accountKey, it.key )
                            .setMountOptions(it.value.mountOptions)

                    mountConfigs << new MountConfiguration().setAzureFileShareConfiguration(shareConfig)
                } else {
                    throw new IllegalArgumentException("Cannot mount a null File Share")
                }
            }
            poolParams.setMountConfiguration(mountConfigs)
        }

        // autoscale
        if( spec.opts.autoScale ) {
            log.debug "Creating autoscale pool with id: ${spec.poolId}; vmCount=${spec.opts.vmCount}; maxVmCount=${spec.opts.maxVmCount}; interval=${spec.opts.scaleInterval}"
            final interval = spec.opts.scaleInterval.seconds as int
            poolParams
                    .setEnableAutoScale(true)
                    .setAutoScaleEvaluationInterval( Duration.of(interval, ChronoUnit.SECONDS) )
                    .setAutoScaleFormula(scaleFormula(spec.opts))
        }
        else if( spec.opts.lowPriority ) {
            log.debug "Creating low-priority pool with id: ${spec.poolId}; vmCount=${spec.opts.vmCount};"
            poolParams
                    .setTargetLowPriorityNodes(spec.opts.vmCount)
        }
        else  {
            log.debug "Creating fixed pool with id: ${spec.poolId}; vmCount=${spec.opts.vmCount};"
            poolParams
                    .setTargetDedicatedNodes(spec.opts.vmCount)
        }

        apply(() -> client.createPool(poolParams))
    }

    protected String scaleFormula(AzPoolOpts opts) {
        final target = opts.lowPriority ? 'TargetLowPriorityNodes' : 'TargetDedicatedNodes'
        // https://docs.microsoft.com/en-us/azure/batch/batch-automatic-scaling
        final DEFAULT_FORMULA = """
            // Get pool lifetime since creation.
            lifespan = time() - time("{{poolCreationTime}}");
            interval = TimeInterval_Minute * {{scaleInterval}};
            
            // Compute the target nodes based on pending tasks.
            // \$PendingTasks == The sum of \$ActiveTasks and \$RunningTasks
            \$samples = \$PendingTasks.GetSamplePercent(interval);
            \$tasks = \$samples < 70 ? max(0, \$PendingTasks.GetSample(1)) : max( \$PendingTasks.GetSample(1), avg(\$PendingTasks.GetSample(interval)));
            \$targetVMs = \$tasks > 0 ? \$tasks : max(0, \$TargetDedicatedNodes/2);
            targetPoolSize = max(0, min(\$targetVMs, {{maxVmCount}}));
            
            // For first interval deploy 1 node, for other intervals scale up/down as per tasks.
            \$${target} = lifespan < interval ? {{vmCount}} : targetPoolSize;
            \$NodeDeallocationOption = taskcompletion;
            """.stripIndent(true)

        final scaleFormula = opts.scaleFormula ?: DEFAULT_FORMULA
        final vars = poolCreationBindings(opts, Instant.now())
        final result = new MustacheTemplateEngine().render(scaleFormula, vars)
        log.debug "Pool autoscale formula:\n$result"
        return result
    }

    protected Map poolCreationBindings(AzPoolOpts opts, Instant time) {
        final vars = new HashMap<String, String>()
        vars.scaleInterval = opts.scaleInterval.minutes
        vars.vmCount = opts.vmCount
        vars.maxVmCount = opts.maxVmCount
        vars.poolCreationTime = time.truncatedTo(ChronoUnit.MICROS).toString()
        return vars
    }

    void deleteTask(AzTaskKey key) {
        apply(() -> client.deleteTask(key.jobId, key.taskId))
    }

    /**
     * Set all jobs to terminate on completion.
     */
    protected void terminateJobs() {
        for( String jobId : allJobIds.values() ) {
            try {
                log.trace "Setting Azure job ${jobId} to terminate on completion"

                final job = apply(() -> client.getJob(jobId))
                final poolInfo = job.poolInfo

                final jobParameter = new BatchJobUpdateContent()
                        .setOnAllTasksComplete(OnAllBatchTasksComplete.TERMINATE_JOB)
                        .setPoolInfo(poolInfo)

                apply(() -> client.updateJob(jobId, jobParameter))
            }
            catch (Exception e) {
                log.warn "Unable to terminate Azure Batch job ${jobId} - Reason: ${e.message ?: e}"
            }
        }
    }

    protected void cleanupJobs() {
        for( String jobId : allJobIds.values() ) {
            try {
                log.trace "Deleting Azure job ${jobId}"
                apply(() -> client.deleteJob(jobId))
            }
            catch (Exception e) {
                log.warn "Unable to delete Azure Batch job ${jobId} - Reason: ${e.message ?: e}"
            }
        }
    }

    protected void cleanupPools() {
        for( String poolId : allPools.keySet() ) {
            try {
                apply(() -> client.deletePool(poolId))
            }
            catch (Exception e) {
                log.warn "Unable to delete Azure Batch pool ${poolId} - Reason: ${e.message ?: e}"
            }
        }
    }

    protected UserIdentity userIdentity(boolean  privileged, String runAs, AutoUserScope scope) {
        UserIdentity identity = new UserIdentity()
        if (runAs) {
            identity.setUsername(runAs)
        } else  {
            identity.setAutoUser(new AutoUserSpecification()
                    .setElevationLevel(privileged ? ElevationLevel.ADMIN : ElevationLevel.NON_ADMIN)
                    .setScope(scope))
        }
        return identity
    }

    @Override
    void close() {
        // terminate all jobs to prevent them from occupying quota
        if( config.batch().terminateJobsOnCompletion ) {
            terminateJobs()
        }

        // delete all jobs
        if( config.batch().deleteJobsOnCompletion ) {
            cleanupJobs()
        }

        // delete all autopools
        if( config.batch().canCreatePool() && config.batch().deletePoolsOnCompletion ) {
            cleanupPools()
        }
    }

    /**
     * Creates a retry policy using the configuration specified by {@link nextflow.cloud.azure.config.AzRetryConfig}
     *
     * @param cond A predicate that determines when a retry should be triggered
     * @return The {@link RetryPolicy} instance
     */
    protected <T> RetryPolicy<T> retryPolicy(Predicate<? extends Throwable> cond) {
        final cfg = config.retryConfig()
        final listener = new EventListener<ExecutionAttemptedEvent<T>>() {
            @Override
            void accept(ExecutionAttemptedEvent<T> event) throws Throwable {
                log.debug("Azure TooManyRequests response error - attempt: ${event.attemptCount}; reason: ${event.lastFailure.message}")
            }
        }
        return RetryPolicy.<T>builder()
                .handleIf(cond)
                .withBackoff(cfg.delay.toMillis(), cfg.maxDelay.toMillis(), ChronoUnit.MILLIS)
                .withMaxAttempts(cfg.maxAttempts)
                .withJitter(cfg.jitter)
                .onRetry(listener)
                .build()
    }

    final private static List<Integer> RETRY_CODES = List.of(408, 429, 500, 502, 503, 504)

    /**
     * Carry out the invocation of the specified action using a retry policy
     * when {@code TooManyRequests} Azure Batch error is returned
     *
     * @param action A {@link CheckedSupplier} instance modeling the action to be performed in a safe manner
     * @return The result of the supplied action
     */
    protected <T> T apply(CheckedSupplier<T> action) {
        // define the retry condition
        final cond = new Predicate<? extends Throwable>() {
            @Override
            boolean test(Throwable t) {
                if( t instanceof HttpResponseException && t.response.statusCode in RETRY_CODES )
                    return true
                if( t instanceof IOException || t.cause instanceof IOException )
                    return true
                if( t instanceof TimeoutException || t.cause instanceof TimeoutException )
                    return true
                return false
            }
        }
        // create the retry policy object
        final policy = retryPolicy(cond)
        // apply the action with
        return Failsafe.with(policy).get(action)
    }
}
