/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.k8s

import java.nio.file.Path
import java.util.concurrent.atomic.AtomicInteger

import groovy.util.logging.Slf4j
import nextflow.container.DockerBuilder
import nextflow.executor.BashWrapperBuilder
import nextflow.k8s.client.K8sClient
import nextflow.k8s.client.K8sResponseException
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.trace.TraceRecord
import nextflow.util.PathTrie
/**
 * Implements the {@link TaskHandler} interface for kubenetes jobs
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class K8sTaskHandler extends TaskHandler {

    static private AtomicInteger VOLUMES = new AtomicInteger()

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

    private Map k8sConfig

    private K8sExecutor executor

    K8sTaskHandler( TaskRun task, K8sExecutor executor ) {
        super(task)
        this.executor = executor
        this.client = executor.client
        this.k8sConfig = executor.k8sConfig
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

    protected boolean getAutoMountHostPaths() {
        def flag = k8sConfig.autoMountHostPaths
        if( flag == null ) {
            // enable auto-mount when kubernetes client is created from
            // a external config file
            flag = !client.config.isFromCluster
        }

        return flag
    }

    protected Map<String,Map> getVolumeClaims() {
        if( !(k8sConfig.volumeClaims instanceof Map) )
            return Collections.emptyMap()

        (Map<String,Map>)k8sConfig.volumeClaims
    }

    protected List<String> getContainerMounts() {

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

        final fixOwnership = builder.fixOwnership()
        final cmd = new ArrayList(new ArrayList(BashWrapperBuilder.BASH)) << TaskRun.CMD_RUN

        final clientConfig = client.config
        final req = [:]
        req.podName = getSyntheticPodName(task)
        req.namespace = clientConfig.namespace
        req.serviceAccount = clientConfig.serviceAccount
        req.imageName = task.container
        req.command = cmd
        req.workDir = task.workDir.toString()
        req.labels = getLabels(task)
        req.env = fixOwnership ? [NXF_OWNER: getOwner()] : Collections.emptyMap()

        // add computing resources
        final taskCfg = task.getConfig()
        final cpus = taskCfg.getCpus()
        final mem = taskCfg.getMemory()
        if( cpus > 1 ) req.cpus = cpus
        if( mem ) req.memory = "${mem.toMega()}Mi"

        // storage
        def claims = getVolumeClaims()
        if( claims ) req.volumeClaims = claims

        final hostMounts = autoMountHostPaths ? getContainerMounts() : Collections.emptyList()
        if( hostMounts ) {
            def mounts = [:]
            hostMounts.each {  path -> mounts[path]=path }
            req.hostMounts = mounts
        }

        K8sHelper.createPodSpec( req )
    }

    protected Map getLabels(TaskRun task) {
        Map result = [:]
        if( executor.getK8sConfig().labels instanceof Map ) {
            executor.getK8sConfig().labels.each { k,v -> result.put(k,sanitize0(v)) }
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

    protected Path yamlDebugPath() {
        boolean debug = executor.getK8sConfig()?.debug?.yaml?.toString() == 'true'
        return debug ? task.workDir.resolve('.command.yaml') : null
    }

    /**
     * @return Retrieve the submitted pod state
     */
    protected Map getState() {
        final now = System.currentTimeMillis()
        final delta =  now - timestamp; timestamp = now
        state && delta < 1_000 ? state : (state = client.podState(podName))
    }

    @Override
    boolean checkIfRunning() {
        if( !podName ) throw new IllegalStateException("Missing K8s pod name -- cannot check if running")
        def state = getState()
        return state.running != null
    }

    @Override
    boolean checkIfCompleted() {
        if( !podName ) throw new IllegalStateException("Missing K8s pod name - cannot check if complete")
        def state = getState()
        if( state.terminated ) {
            // finalize the task
            task.exitStatus = readExitFile()
            task.stdout = outputFile
            task.stderr = errorFile
            status = TaskStatus.COMPLETED
            deletePodIfSuccessful(task)
            return true
        }

        return false
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
        if( podName ) {
            log.trace "[K8s] deleting pod name=$podName"
            client.podDelete(podName)
        }
        else {
            log.debug "[K8s] Oops.. invalid delete action"
        }
    }

    protected void deletePodIfSuccessful(TaskRun task) {
        if( !podName )
            return

        if( executor.getK8sConfig().cleanup == false )
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
