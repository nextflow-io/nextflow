/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.container.DockerBuilder
import nextflow.exception.ProcessSubmitException
import nextflow.executor.BashWrapperBuilder
import nextflow.k8s.client.K8sClient
import nextflow.k8s.client.K8sResponseException
import nextflow.k8s.model.PodEnv
import nextflow.k8s.model.PodOptions
import nextflow.k8s.model.PodSpecBuilder
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.trace.TraceRecord
import nextflow.util.PathTrie
/**
 * Implements the {@link TaskHandler} interface for Kubernetes jobs
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


    private K8sClient client

    private String podName

    private K8sWrapperBuilder builder

    private Path outputFile

    private Path errorFile

    private Path exitFile

    private Map state

    private long timestamp

    private K8sExecutor executor

    K8sTaskHandler( TaskRun task, K8sExecutor executor ) {
        super(task)
        this.executor = executor
        this.client = executor.client
        this.outputFile = task.workDir.resolve(TaskRun.CMD_OUTFILE)
        this.errorFile = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        this.exitFile = task.workDir.resolve(TaskRun.CMD_EXIT)
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

    protected K8sConfig getK8sConfig() { executor.getK8sConfig() }

    protected List<String> getContainerMounts() {

        if( !k8sConfig.getAutoMountHostPaths() ) {
            return Collections.<String>emptyList()
        }

        // get input files paths
        final paths = DockerBuilder.inputFilesToPaths(builder.getResolvedInputs())
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
            throw  new ProcessSubmitException("Failed to submit K8s job -- Cause: ${e.message ?: e}", e)
        }
    }

    protected Map newSubmitRequest0(TaskRun task, String imageName) {

        final fixOwnership = builder.fixOwnership()
        final cmd = new ArrayList(new ArrayList(BashWrapperBuilder.BASH)) << TaskRun.CMD_RUN
        final taskCfg = task.getConfig()

        final clientConfig = client.config
        final builder = new PodSpecBuilder()
            .withImageName(imageName)
            .withPodName(getSyntheticPodName(task))
            .withCommand(cmd)
            .withWorkDir(task.workDir)
            .withNamespace(clientConfig.namespace)
            .withServiceAccount(clientConfig.serviceAccount)
            .withLabels(getLabels(task))
            .withPodOptions(getPodOptions())

        // note: task environment is managed by the task bash wrapper
        // do not add here -- see also #680
        if( fixOwnership )
            builder.withEnv(PodEnv.value('NXF_OWNER', getOwner()))

        // add computing resources
        final cpus = taskCfg.get('cpus') as Integer
        final mem = taskCfg.getMemory()
        if( cpus )
            builder.withCpus(cpus as int)
        if( mem )
            builder.withMemory(mem)

        final List<String> hostMounts = getContainerMounts()
        for( String mount : hostMounts ) {
            builder.withHostMount(mount,mount)
        }

        return builder.build()
    }

    protected PodOptions getPodOptions() {
        // merge the pod options provided in the k8s config
        // with the ones in process config
        def opt1 = k8sConfig.getPodOptions()
        def opt2 = task.getConfig().getPodOptions()
        return opt1 + opt2
    }


    protected Map getLabels(TaskRun task) {
        Map result = [:]
        def labels = k8sConfig.getLabels()
        if( labels ) {
            labels.each { k,v -> result.put(k,sanitize0(v)) }
        }
        result.app = 'nextflow'
        result.runName = sanitize0(getRunName())
        result.taskName = sanitize0(task.getName())
        result.processName = sanitize0(task.getProcessor().getName())
        result.sessionId = sanitize0("uuid-${executor.getSession().uniqueId}")
        return result
    }

    /**
     * Valid label must be an empty string or consist of alphanumeric characters, '-', '_' or '.',
     * and must start and end with an alphanumeric character.
     *
     * @param value
     * @return
     */

    protected String sanitize0( value ) {
        def str = String.valueOf(value)
        str = str.replaceAll(/[^a-zA-Z0-9\.\_\-]+/, '_')
        str = str.replaceAll(/^[^a-zA-Z]+/, '')
        str = str.replaceAll(/[^a-zA-Z0-9]+$/, '')
        return str
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
        final resp = client.podCreate(req, yamlDebugPath())

        if( !resp.metadata?.name )
            throw new K8sResponseException("Missing created pod name", resp)
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
        final delta =  now - timestamp;
        if( !state || delta >= 1_000) {
            def newState = client.podState(podName)
            if( newState ) {
                state = newState
                timestamp = now
            }
        }
        return state
    }

    @Override
    boolean checkIfRunning() {
        if( !podName ) throw new IllegalStateException("Missing K8s pod name -- cannot check if running")
        def state = getState()
        return state && state.running != null
    }

    @Override
    boolean checkIfCompleted() {
        if( !podName ) throw new IllegalStateException("Missing K8s pod name - cannot check if complete")
        def state = getState()
        if( state && state.terminated ) {
            // finalize the task
            task.exitStatus = readExitFile()
            task.stdout = outputFile
            task.stderr = errorFile
            status = TaskStatus.COMPLETED
            savePodLogOnError(task)
            deletePodIfSuccessful(task)
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
            final stream = client.podLog(podName)
            Files.copy(stream, task.workDir.resolve(TaskRun.CMD_LOG))
        }
        catch( Exception e ) {
            log.warn "Failed to copy log for pod $podName", e
        }
    }

    protected int readExitFile() {
        try {
            exitFile.text as Integer
        }
        catch( Exception e ) {
            log.debug "[K8s] Cannot read exitstatus for task: `$task.name`", e
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
            log.trace "[K8s] deleting pod name=$podName"
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
            client.podDelete(podName)
        }
        catch( Exception e ) {
            log.warn "Unable to cleanup pod: $podName -- see the log file for details", e
        }
    }


    TraceRecord getTraceRecord() {
        final result = super.getTraceRecord()
        result.put('native_id', podName)
        return result
    }

}
