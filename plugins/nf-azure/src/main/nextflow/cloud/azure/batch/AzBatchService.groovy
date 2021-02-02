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

import com.microsoft.azure.batch.protocol.models.StartTask

import java.math.RoundingMode
import java.nio.file.Path
import java.time.Instant

import com.microsoft.azure.batch.BatchClient
import com.microsoft.azure.batch.auth.BatchSharedKeyCredentials
import com.microsoft.azure.batch.protocol.models.BatchErrorException
import com.microsoft.azure.batch.protocol.models.CloudPool
import com.microsoft.azure.batch.protocol.models.CloudTask
import com.microsoft.azure.batch.protocol.models.ContainerConfiguration
import com.microsoft.azure.batch.protocol.models.ImageInformation
import com.microsoft.azure.batch.protocol.models.OutputFile
import com.microsoft.azure.batch.protocol.models.OutputFileBlobContainerDestination
import com.microsoft.azure.batch.protocol.models.OutputFileDestination
import com.microsoft.azure.batch.protocol.models.OutputFileUploadCondition
import com.microsoft.azure.batch.protocol.models.OutputFileUploadOptions
import com.microsoft.azure.batch.protocol.models.PoolAddParameter
import com.microsoft.azure.batch.protocol.models.PoolInformation
import com.microsoft.azure.batch.protocol.models.PoolState
import com.microsoft.azure.batch.protocol.models.ResourceFile
import com.microsoft.azure.batch.protocol.models.TaskAddParameter
import com.microsoft.azure.batch.protocol.models.TaskContainerSettings
import com.microsoft.azure.batch.protocol.models.VirtualMachineConfiguration
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
import nextflow.cloud.azure.config.AzConfig
import nextflow.cloud.azure.config.AzPoolOpts
import nextflow.cloud.azure.nio.AzPath
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.util.CacheHelper
import nextflow.util.MemoryUnit
import nextflow.util.MustacheTemplateEngine
import nextflow.util.Rnd
import org.joda.time.Period
/**
 * Implements Azure Batch operations for Nextflow executor
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AzBatchService implements Closeable {

    static private final long _1GB = 1 << 30

    static final private Map<String,AzVmPoolSpec> allPools = new HashMap<>(50)

    AzConfig config

    Map<TaskProcessor,String> allJobIds = new HashMap<>(50)

    AzBatchService(AzBatchExecutor executor) {
        assert executor
        this.config = executor.config
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

    AzVmType findBestVm(String location, int cpus, MemoryUnit mem, String family) {
        def all = listAllVms(location)
        def scores = new TreeMap<Double,String>()
        for( Map entry : all ) {
            if( !matchType(family, entry.name as String) )
                continue
            def score = computeScore(cpus, mem, entry)
            if( score != null )
                scores.put(score, entry.name as String)
        }
        return getVmType(location, scores.firstEntry().value)
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
        double vmMemGb = (entry.memoryInMb as int) /1024

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

        return Math.sqrt(score);
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

    protected BatchClient createBatchClient() {
        log.debug "[AZURE BATCH] Executor options=${config.batch()}"
        // Create batch client
        if( !config.batch().endpoint )
            throw new IllegalArgumentException("Missing Azure Batch endpoint -- Specify it in the nextflow.config file using the setting 'azure.batch.endpoint'")
        if( !config.batch().accountName )
            throw new IllegalArgumentException("Missing Azure Batch account name -- Specify it in the nextflow.config file using the setting 'azure.batch.accountName'")
        if( !config.batch().accountKey )
            throw new IllegalArgumentException("Missing Azure Batch account key -- Specify it in the nextflow.config file using the setting 'azure.batch.accountKet'")

        final cred = new BatchSharedKeyCredentials(config.batch().endpoint, config.batch().accountName, config.batch().accountKey)
        final client = BatchClient.open(cred)
        final sess = Global.session as Session
        sess.onShutdown { client.protocolLayer().restClient().close() }
        return client
    }

    AzTaskKey submitTask(TaskRun task) {
        final poolId = getOrCreatePool(task)
        final jobId = getOrCreateJob(poolId, task)
        runTask(poolId, jobId, task)
    }

    CloudTask getTask(AzTaskKey key) {
        return client.taskOperations().getTask(key.jobId, key.taskId)
    }

    void terminate(AzTaskKey key) {
        client.taskOperations().terminateTask(key.jobId, key.taskId)
    }

    synchronized String getOrCreateJob(String poolId, TaskRun task) {
        final mapKey = task.processor
        if( allJobIds.containsKey(mapKey)) {
            return allJobIds[mapKey]
        }
        // create a batch job
        final jobId = makeJobId(task)
        final poolInfo = new PoolInformation()
                            .withPoolId(poolId)
        client
            .jobOperations()
            .createJob(jobId, poolInfo)
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

        final key = "job-${Rnd.hex()}-${name}"
        return key.size()>64 ? key.substring(0,64) : key
    }

    AzTaskKey runTask(String poolId, String jobId, TaskRun task) {
        assert poolId, 'Missing Azure Batch poolId argument'
        assert jobId, 'Missing Azure Batch jobId argument'
        assert task, 'Missing Azure Batch task argument'

        final sas = config.storage().sasToken
        if( !sas )
            throw new IllegalArgumentException("Missing Azure Blob storage SAS token")

        final container = task.config.container as String
        if( !container )
            throw new IllegalArgumentException("Missing container image for process: $task.name")

        final taskId = "nf-${task.hash.toString()}"
        final containerOpts = new TaskContainerSettings()
                .withImageName(container)
                // mount host certificates otherwise `azcopy fails
                .withContainerRunOptions('-v /etc/ssl/certs:/etc/ssl/certs:ro -v /etc/pki:/etc/pki:ro')

        final pool = allPools.get(poolId)
        if( !pool )
            throw new IllegalStateException("Missing Azure Batch pool spec with id: $poolId")

        final slots = computeSlots(task, pool)
        log.debug "[AZURE BATCH] Submitting task: $taskId, cpus=${task.config.getCpus()}, mem=${task.config.getMemory()?:'-'}, slots: $slots"

        final taskToAdd = new TaskAddParameter()
                .withId(taskId)
                .withContainerSettings(containerOpts)
                .withCommandLine("sh -c 'bash ${TaskRun.CMD_RUN} 2>&1 | tee ${TaskRun.CMD_LOG}'")
                .withResourceFiles(resourceFileUrls(task,sas))
                .withOutputFiles(outputFileUrls(task, sas))
                .withRequiredSlots(slots)

        client.taskOperations().createTask(jobId, taskToAdd)
        return new AzTaskKey(jobId, taskId)
    }

    protected List<ResourceFile> resourceFileUrls(TaskRun task, String sas) {
        final cmdRun = (AzPath) task.workDir.resolve(TaskRun.CMD_RUN)
        final cmdScript = (AzPath) task.workDir.resolve(TaskRun.CMD_SCRIPT)

        final resFiles = new ArrayList(10)

        resFiles << new ResourceFile()
                .withHttpUrl(AzHelper.toHttpUrl(cmdRun, sas))
                .withFilePath(TaskRun.CMD_RUN)

        resFiles << new ResourceFile()
                .withHttpUrl(AzHelper.toHttpUrl(cmdScript, sas))
                .withFilePath(TaskRun.CMD_SCRIPT)

        if( task.stdin ) {
            resFiles << new ResourceFile()
                    .withHttpUrl(AzHelper.toHttpUrl(cmdScript, sas))
                    .withFilePath(TaskRun.CMD_INFILE)
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
        final dest = new OutputFileBlobContainerDestination()
                .withContainerUrl(AzHelper.toContainerUrl(targetDir,sas))
                .withPath(target.subpath(1,target.nameCount).toString())

        return new OutputFile()
                .withFilePattern(localPath)
                .withDestination( new OutputFileDestination().withContainer(dest) )
                .withUploadOptions( new OutputFileUploadOptions().withUploadCondition( OutputFileUploadCondition.TASK_COMPLETION) )
    }

    protected ImageInformation getImage(AzPoolOpts opts) {
        List<ImageInformation> images = client.accountOperations().listSupportedImages()

        for (ImageInformation it : images) {
            if( it.osType() != opts.osType )
                continue
            if( it.verificationType() != opts.verification )
                continue
            if( !it.imageReference().publisher().equalsIgnoreCase(opts.publisher) )
                continue
            if( it.imageReference().offer().equalsIgnoreCase(opts.offer) )
                return it
        }

        throw new IllegalStateException("Cannot find a matching VM image with publister=$opts.publisher; offer=$opts.offer; OS type=$opts.osType; verification type=$opts.verification")
    }

    protected AzVmPoolSpec specFromConfigPool(String poolId) {

        def opts = config.batch().pool(poolId)
        if( !opts )
            throw new IllegalArgumentException("Cannot find Azure Batch config for pool: $poolId")

        def type = getVmType(config.batch().location, opts.vmType)
        if( !type )
            throw new IllegalArgumentException("Cannot find Azure Batch VM type '$poolId' - Check pool definition $poolId in the Nextflow config file")

        new AzVmPoolSpec(poolId: poolId, vmType: type, opts: opts)
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
            throw new IllegalArgumentException("Missing Azure Batch VM type")

        final vmType = guessBestVm(loc, cpus, mem, type)
        if( !vmType ) {
            def msg = "Cannot find a VM matching this requirements: type=$type, cpus=${cpus}, mem=${mem?:'-'}, location=${loc}"
            throw new IllegalArgumentException(msg)
        }

        final key = CacheHelper.hasher([vmType.name, opts]).hash().toString()
        final poolId = "nf-pool-$key-$vmType.name"
        return new AzVmPoolSpec(poolId: poolId, vmType: vmType, opts: opts)
    }

    protected void checkPool(CloudPool pool, AzVmPoolSpec spec) {
        if( pool.state() != PoolState.ACTIVE ) {
            throw new IllegalStateException("Azure Batch pool '${pool.id()}' not in active state")
        }
        else if (pool.resizeErrors() && pool.currentDedicatedNodes()==0 ) {
            throw new IllegalStateException("Azure Batch pool '${pool.id()}' has resize errors")
        }
        if( pool.taskSlotsPerNode() != spec.vmType.numberOfCores ) {
            throw new IllegalStateException("Azure Batch pool '${pool.id()}' slots per node does not match the VM num cores (slots: ${pool.taskSlotsPerNode()}, cores: ${spec.vmType.numberOfCores})")
        }
    }

    protected void checkPoolId(String poolId) {
        if( !poolId.matches(/^[\w\-]+$/) )
            throw new IllegalArgumentException("Invalid Azure Batch pool Id '$poolId' - It can only contains alphanumeric, hyphen and undershore characters")
    }

    protected AzVmPoolSpec specForTask(TaskRun task) {
        String poolId = null
        if( !config.batch().autoPoolMode ) {
            // the process queue is used as poolId
            poolId = task.config.queue as String
            if( !poolId ) {
                throw new IllegalArgumentException("No Azure Batch pool was specified for task '${task.name}' - Either specify the pool name using the 'queue' diretive or enable the 'autoPoolMode' option")
            }
            // sanity check
            checkPoolId(poolId)
            // check if cached
            if( allPools.containsKey(poolId) ) {
                return allPools.get(poolId)
            }
        }

        return poolId
                ? specFromConfigPool(poolId)
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
                throw new IllegalArgumentException("Can't find Azure Batch pool '$spec.poolId' - Make sure it exists or enablethe use `allowPoolCreation=true` in the nextflow config file")
            }
        }
        else {
            checkPool(pool, spec)
        }

        // add to the list of pool ids
        allPools[spec.poolId] = spec

        return spec.poolId
    }


    protected CloudPool getPool(String poolId) {
        try {
            return client.poolOperations().getPool(poolId)
        }
        catch (BatchErrorException e) {
            if( e.response().code() == 404 ) {
                // not found
                return null
            }
            throw e
        }
    }

    protected VirtualMachineConfiguration poolVmConfig(AzPoolOpts opts) {
        /**
         * A container configuration must be provided for a task to run in a specific container.
         * Such container can be pre-fetched on VM creation or when running the task
         *
         * https://github.com/MicrosoftDocs/azure-docs/blob/master/articles/batch/batch-docker-container-workloads.md#:~:text=Run%20container%20applications%20on%20Azure,compatible%20containers%20on%20the%20nodes.
         */
        final containerConfig = new ContainerConfiguration();
        final image = getImage(opts)

        new VirtualMachineConfiguration()
                .withNodeAgentSKUId(image.nodeAgentSKUId())
                .withImageReference(image.imageReference())
                .withContainerConfiguration(containerConfig)
    }

    protected void createPool(AzVmPoolSpec spec) {

        def resourceFiles = new ArrayList(10)
        
        resourceFiles << new ResourceFile()
                .withHttpUrl("https://nf-xpack.s3-eu-west-1.amazonaws.com/azcopy/linux_amd64_10.8.0/azcopy")
                .withFilePath('azcopy')

        def poolStartTask = new StartTask()
                .withCommandLine('bash -c "chmod +x azcopy && mkdir \$AZ_BATCH_NODE_SHARED_DIR/bin/ && cp azcopy \$AZ_BATCH_NODE_SHARED_DIR/bin/" ')
                .withResourceFiles(resourceFiles)


        final poolParams = new PoolAddParameter()
                .withId(spec.poolId)
                .withVirtualMachineConfiguration(poolVmConfig(spec.opts))
                // https://docs.microsoft.com/en-us/azure/batch/batch-pool-vm-sizes
                .withVmSize(spec.vmType.name)
                // same as the num ofd cores
                // https://docs.microsoft.com/en-us/azure/batch/batch-parallel-node-tasks
                .withTaskSlotsPerNode(spec.vmType.numberOfCores)
                .withStartTask(poolStartTask)

        if( spec.opts.autoScale ) {
            poolParams
                    .withEnableAutoScale(true)
                    .withAutoScaleEvaluationInterval( new Period().withSeconds(300) ) // cannot be smaller
                    .withAutoScaleFormula(scaleFormula(spec.opts))
        }
        else {
            poolParams
                    .withTargetDedicatedNodes(spec.opts.vmCount)
        }

        client.poolOperations().createPool(poolParams)
    }

    protected String scaleFormula(AzPoolOpts opts) {
        // https://docs.microsoft.com/en-us/azure/batch/batch-automatic-scaling
        def DEFAULT_FORMULA = '''\
                $TargetDedicatedNodes = {{vmCount}};
                lifespan         = time() - time("{{now}}");
                span             = TimeInterval_Minute * 60;
                startup          = TimeInterval_Minute * 10;
                ratio            = 50;
                $TargetDedicatedNodes = (lifespan > startup ? (max($RunningTasks.GetSample(span, ratio), $ActiveTasks.GetSample(span, ratio)) == 0 ? 0 : $TargetDedicatedNodes) : {{vmCount}});
                '''.stripIndent()

        final scaleFormula = opts.scaleFormula ?: DEFAULT_FORMULA
        final vars = new HashMap<String,String>()
        vars.vmCount = opts.vmCount
        vars.now = Instant.now().toString()
        return new MustacheTemplateEngine().render(scaleFormula, vars)
    }

    void deleteTask(AzTaskKey key) {
        client.taskOperations().deleteTask(key.jobId, key.taskId)
    }

    protected void cleanupJobs() {
        for( Map.Entry<TaskProcessor,String> entry : allJobIds ) {
            final proc = entry.key
            final jobId = entry.value
            if( proc.hasErrors() ) {
                log.debug "Preserving Azure job with error: ${jobId}"
                continue
            }

            try {
                log.trace "Deleting Azure job ${jobId}"
                client.jobOperations().deleteJob(jobId)
            }
            catch (Exception e) {
                log.warn "Unable to delete Azure Batch job ${jobId} - Reason: ${e.message ?: e}"
            }
        }
    }

    protected void cleanupPools() {
        for( String poolId : allPools.keySet()) {
            try {
                client.poolOperations().deletePool(poolId)
            }
            catch (Exception e) {
                log.warn "Unable to delete Azure Batch pool ${poolId} - Reason: ${e.message ?: e}"
            }
        }
    }

    @Override
    void close() {
        // cleanup app successful jobs
        if( config.batch().deleteJobsOnCompletion!=Boolean.FALSE ) {
            cleanupJobs()
        }
        if( config.batch().canCreatePool() && config.batch().deletePoolsOnCompletion ) {
            cleanupPools()
        }
    }
}
