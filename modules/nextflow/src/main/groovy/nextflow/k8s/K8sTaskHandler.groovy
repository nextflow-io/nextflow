/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.k8s

import java.nio.file.Files
import java.nio.file.Path
import java.time.Instant
import java.time.format.DateTimeFormatter

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.container.DockerBuilder
import nextflow.exception.NodeTerminationException
import nextflow.exception.ProcessSubmitException
import nextflow.executor.BashWrapperBuilder
import nextflow.k8s.client.K8sClient
import nextflow.k8s.client.K8sResponseException
import nextflow.k8s.model.PodEnv
import nextflow.k8s.model.PodOptions
import nextflow.k8s.model.PodSpecBuilder
import nextflow.k8s.model.ResourceType
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.trace.TraceRecord
import nextflow.util.Escape
import nextflow.util.PathTrie
/**
 * Implements the {@link TaskHandler} interface for Kubernetes pods
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class K8sTaskHandler extends TaskHandler {

    @Lazy
    static private final String OWNER = {
        if( System.getenv('NXF_OWNER') ) {
            return System.getenv('NXF_OWNER')
        }
        else {
            def p = ['bash','-c','echo -n $(id -u):$(id -g)'].execute();
            p.waitFor()
            return p.text
        }

    } ()

    private ResourceType resourceType = ResourceType.Pod

    private K8sClient client

    private String podName

    private K8sWrapperBuilder builder

    private Path outputFile

    private Path errorFile

    private Path exitFile

    private Map state

    private long timestamp

    private K8sExecutor executor

    private String runsOnNode = null

    K8sTaskHandler( TaskRun task, K8sExecutor executor ) {
        super(task)
        this.executor = executor
        this.client = executor.client
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
        this.resourceType = executor.k8sConfig.useJobResource() ? ResourceType.Job : ResourceType.Pod
    }

    /** only for testing -- do not use */
    protected K8sTaskHandler() {

    }

    /**
     * @return The workflow execution unique run name
     */
    protected String getRunName() {
        executor.session.runName
    }

    protected String getPodName() {
        return podName
    }

    protected K8sConfig getK8sConfig() { executor.getK8sConfig() }

    protected boolean useJobResource() { resourceType==ResourceType.Job }

    protected List<String> getContainerMounts() {

        if( !k8sConfig.getAutoMountHostPaths() ) {
            return Collections.<String>emptyList()
        }

        // get input files paths
        final paths = DockerBuilder.inputFilesToPaths(builder.getInputFiles())
        final binDir = builder.binDir
        final workDir = builder.workDir
        // add standard paths
        if( binDir ) paths << binDir
        if( workDir ) paths << workDir

        def trie = new PathTrie()
        paths.each { trie.add(it) }

        // defines the mounts
        trie.longest()
    }

    protected K8sWrapperBuilder createBashWrapper(TaskRun task) {
        new K8sWrapperBuilder(task)
    }

    protected String getSyntheticPodName(TaskRun task) {
        "nf-${task.hash}"
    }

    protected String getOwner() { OWNER }

    /**
     * Creates a Pod specification that executed that specified task
     *
     * @param task A {@link TaskRun} instance representing the task to execute
     * @return A {@link Map} object modeling a Pod specification
     */

    protected Map newSubmitRequest(TaskRun task) {
        def imageName = task.container
        if( !imageName )
            throw new ProcessSubmitException("Missing container image for process `$task.processor.name`")

        try {
            newSubmitRequest0(task, imageName)
        }
        catch( Throwable e ) {
            throw  new ProcessSubmitException("Failed to submit K8s ${resourceType.lower()} -- Cause: ${e.message ?: e}", e)
        }
    }

    protected boolean preserveContainerEntrypoint() {
        return executor.getK8sConfig().preserveContainerEntrypoint()
    }

    protected Map newSubmitRequest0(TaskRun task, String imageName) {

        final fixOwnership = builder.fixOwnership()
        final launcher = new ArrayList(new ArrayList(BashWrapperBuilder.BASH)) << "${Escape.path(task.workDir)}/${TaskRun.CMD_RUN}".toString()
        final taskCfg = task.getConfig()

        final clientConfig = client.config
        final builder = new PodSpecBuilder()
            .withImageName(imageName)
            .withPodName(getSyntheticPodName(task))
            .withNamespace(clientConfig.namespace)
            .withServiceAccount(clientConfig.serviceAccount)
            .withLabels(getLabels(task))
            .withAnnotations(getAnnotations())
            .withPodOptions(getPodOptions())

        // when `preserveEntrypoint` is true the launcher is run via `args` instead of `command`
        // to not override the container entrypoint
        if( preserveContainerEntrypoint() ) {
            builder.withArgs(launcher)
        }
        else {
            builder.withCommand(launcher)
        }

        // note: task environment is managed by the task bash wrapper
        // do not add here -- see also #680
        if( fixOwnership )
            builder.withEnv(PodEnv.value('NXF_OWNER', getOwner()))

        // add computing resources
        final cpus = taskCfg.getCpuUnits()
        final mem = taskCfg.getMemory()
        final acc = taskCfg.getAccelerator()
        if( cpus )
            builder.withCpus(cpus)
        if( mem )
            builder.withMemory(mem)
        if( acc )
            builder.withAccelerator(acc)

        final List<String> hostMounts = getContainerMounts()
        for( String mount : hostMounts ) {
            builder.withHostMount(mount,mount)
        }

        if ( taskCfg.time ) {
            final duration = taskCfg.getTime()
            builder.withActiveDeadline(duration.toSeconds() as int)
        }

        return useJobResource()
            ? builder.buildAsJob()
            : builder.build()
    }

    protected PodOptions getPodOptions() {
        // merge the pod options provided in the k8s config
        // with the ones in process config
        def opt1 = k8sConfig.getPodOptions()
        def opt2 = task.getConfig().getPodOptions()
        return opt1 + opt2
    }


    protected Map<String,String> getLabels(TaskRun task) {
        def result = new LinkedHashMap<String,String>(10)
        def labels = k8sConfig.getLabels()
        if( labels ) {
            result.putAll(labels)
        }
        result.app = 'nextflow'
        result.runName = getRunName()
        result.taskName = task.getName()
        result.processName = task.getProcessor().getName()
        result.sessionId = "uuid-${executor.getSession().uniqueId}" as String
        return result
    }

    protected Map getAnnotations() {
        k8sConfig.getAnnotations()
    }

    /**
     * Creates a new K8s pod executing the associated task
     */
    @Override
    @CompileDynamic
    void submit() {
        builder = createBashWrapper(task)
        builder.build()

        final req = newSubmitRequest(task)
        final resp = useJobResource()
                ? client.jobCreate(req, yamlDebugPath())
                : client.podCreate(req, yamlDebugPath())

        if( !resp.metadata?.name )
            throw new K8sResponseException("Missing created ${resourceType.lower()} name", resp)
        this.podName = resp.metadata.name
        this.status = TaskStatus.SUBMITTED
    }

    @CompileDynamic
    protected Path yamlDebugPath() {
        boolean debug = k8sConfig.getDebug().getYaml()
        return debug ? task.workDir.resolve('.command.yaml') : null
    }

    /**
     * @return Retrieve the submitted pod state
     */
    protected Map getState() {
        final now = System.currentTimeMillis()
        try {
            final delta =  now - timestamp;
            if( !state || delta >= 1_000) {
                def newState = useJobResource()
                        ? client.jobState(podName)
                        : client.podState(podName)
                if( newState ) {
                   log.trace "[K8s] Get ${resourceType.lower()}=$podName state=$newState"
                   state = newState
                   timestamp = now
                }
            }
            return state
        } 
        catch (NodeTerminationException e) {
            // create a synthetic `state` object adding an extra `nodeTermination`
            // attribute to return the NodeTerminationException error to the caller method
            final instant = Instant.now()
            final result = new HashMap(10)
            result.terminated = [startedAt:instant.toString(), finishedAt:instant.toString()]
            result.nodeTermination = e
            timestamp = now
            state = result
            return state
        }
    }

    @Override
    boolean checkIfRunning() {
        if( !podName ) throw new IllegalStateException("Missing K8s ${resourceType.lower()} name -- cannot check if running")
        if(isSubmitted()) {
            def state = getState()
            // include `terminated` state to allow the handler status to progress
            if (state && (state.running != null || state.terminated)) {
                status = TaskStatus.RUNNING
                determineNode()
                return true
            }
        }
        return false
    }

    long getEpochMilli(String timeString) {
        final time = DateTimeFormatter.ISO_INSTANT.parse(timeString)
        return Instant.from(time).toEpochMilli()
    }

    /**
     * Update task start and end times based on pod timestamps.
     * We update timestamps because it's possible for a task to run  so quickly
     * (less than 1 second) that it skips right over the RUNNING status.
     * If this happens, the startTimeMillis never gets set and remains equal to 0.
     * To make sure startTimeMillis is non-zero we update it with the pod start time.
     * We update completeTimeMillis from the same pod info to be consistent.
     */
    void updateTimestamps(Map terminated) {
        try {
            startTimeMillis = getEpochMilli(terminated.startedAt as String)
            completeTimeMillis = getEpochMilli(terminated.finishedAt as String)
        } catch( Exception e ) {
            log.debug "Failed updating timestamps '${terminated.toString()}'", e
            // Only update if startTimeMillis hasn't already been set.
            // If startTimeMillis _has_ been set, then both startTimeMillis
            // and completeTimeMillis will have been set with the normal
            // TaskHandler mechanism, so there's no need to reset them here.
            if (!startTimeMillis) {
                startTimeMillis = System.currentTimeMillis()
                completeTimeMillis = System.currentTimeMillis()
            }
        }
    }

    @Override
    boolean checkIfCompleted() {
        if( !podName ) throw new IllegalStateException("Missing K8s ${resourceType.lower()} name - cannot check if complete")
        def state = getState()
        if( state && state.terminated ) {
            if( state.nodeTermination instanceof NodeTerminationException ) {
                // kee track of the node termination error
                task.error = (NodeTerminationException) state.nodeTermination
                // mark the task as ABORTED since thr failure is caused by a node failure
                task.aborted = true
            }
            else {
                // finalize the task
                task.exitStatus = readExitFile()
                task.stdout = outputFile
                task.stderr = errorFile
            }
            status = TaskStatus.COMPLETED
            savePodLogOnError(task)
            deletePodIfSuccessful(task)
            updateTimestamps(state.terminated as Map)
            determineNode()
            return true
        }

        return false
    }

    protected void savePodLogOnError(TaskRun task) {
        if( task.isSuccess() )
            return

        if( errorFile && !errorFile.empty() )
            return

        final session = executor.getSession()
        if( session.isAborted() || session.isCancelled() || session.isTerminated() )
            return

        try {
            final stream = useJobResource()
                    ? client.jobLog(podName)
                    : client.podLog(podName)
            Files.copy(stream, task.workDir.resolve(TaskRun.CMD_LOG))
        }
        catch( Exception e ) {
            log.warn "Failed to copy log for ${resourceType.lower()} $podName", e
        }
    }

    protected int readExitFile() {
        try {
            exitFile.text as Integer
        }
        catch( Exception e ) {
            log.debug "[K8s] Cannot read exitstatus for task: `$task.name` | ${e.message}"
            return Integer.MAX_VALUE
        }
    }

    /**
     * Terminates the current task execution
     */
    @Override
    void kill() {
        if( cleanupDisabled() )
            return
        
        if( podName ) {
            log.trace "[K8s] deleting ${resourceType.lower()} name=$podName"
            if ( useJobResource() )
                client.jobDelete(podName)
            else
                client.podDelete(podName)
        }
        else {
            log.debug "[K8s] Oops.. invalid delete action"
        }
    }

    protected boolean cleanupDisabled() {
        !k8sConfig.getCleanup()
    }

    protected void deletePodIfSuccessful(TaskRun task) {
        if( !podName )
            return

        if( cleanupDisabled() )
            return

        if( !task.isSuccess() ) {
            // do not delete successfully executed pods for debugging purpose
            return
        }

        try {
            if ( useJobResource() )
                client.jobDelete(podName)
            else
                client.podDelete(podName)
        }
        catch( Exception e ) {
            log.warn "Unable to cleanup ${resourceType.lower()}: $podName -- see the log file for details", e
        }
    }

    private void determineNode(){
        try {
            if ( k8sConfig.fetchNodeName() && !runsOnNode )
                runsOnNode = client.getNodeOfPod( podName )
        } catch ( Exception e ){
            log.warn ("Unable to fetch pod: $podName its node -- see the log file for details", e)
        }
    }

    TraceRecord getTraceRecord() {
        final result = super.getTraceRecord()
        result.put('native_id', podName)
        result.put( 'hostname', runsOnNode )
        return result
    }

}
