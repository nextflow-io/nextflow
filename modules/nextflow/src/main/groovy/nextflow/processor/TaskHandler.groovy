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

package nextflow.processor

import static nextflow.processor.TaskStatus.*
import static nextflow.trace.TraceRecord.toLong
import static nextflow.trace.TraceRecord.toFloat

import java.nio.file.NoSuchFileException
import java.util.concurrent.atomic.AtomicBoolean

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.trace.TraceRecord
import nextflow.util.TestOnly
/**
 * Actions to handle the underlying job running the user task.
 *
 * <p>
 * Note this types the job in the execution facility (i.e. grid, cloud, etc)
 * NOT the *logical* user task
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class TaskHandler {

    private AtomicBoolean killed = new AtomicBoolean()

    protected TaskHandler(TaskRun task) {
        this.task = task
    }

    @TestOnly
    protected TaskHandler() {}

    /**
     * The task managed by this handler
     */
    TaskRun task

    /**
     * The task managed by this handler
     */
    TaskRun getTask() { task }

    /**
     * Whenever this handle reference a job array task child
     */
    boolean isArrayChild

    TaskHandler withArrayChild(boolean child) {
        this.isArrayChild = child
        return this
    }

    /**
     * Task current status
     */
    volatile TaskStatus status = NEW

    long submitTimeMillis

    long startTimeMillis

    long completeTimeMillis


    /**
     * Model the start transition from {@code #SUBMITTED} to {@code STARTED}
     */
    abstract boolean checkIfRunning()

    /**
     *  Model the start transition from {@code #STARTED} to {@code TERMINATED}
     */
    abstract boolean checkIfCompleted()

    /**
     * Template method implementing the termination of a task execution.
     * This is not mean to be invoked directly. See also {@link #kill()}
     */
    abstract protected void killTask()

    /**
     * Kill a job execution.
     *
     * @see #killTask()
     */
    void kill() {
        if (!killed.getAndSet(true)) {
            killTask()
        }
    }

    /**
     * Submit the task for execution.
     *
     * Note: the underlying execution platform may schedule it in its own queue
     */
    abstract void submit()

    /**
     * Prepare the launcher script.
     *
     * This method is optional. If it is not implemented, the launcher script should
     * be prepared in the submit() method.
     */
    void prepareLauncher() {}

    /**
     * Task status attribute setter.
     *
     * @param status The sask status as defined by {@link TaskStatus}
     */
    void setStatus( TaskStatus status ) {

        // skip if the status is the same aam
        if ( this.status == status || status == null )
            return

        // change the status
        this.status = status
        switch( status ) {
            case SUBMITTED: submitTimeMillis = System.currentTimeMillis(); break
            case RUNNING: startTimeMillis = System.currentTimeMillis(); break
            case COMPLETED: completeTimeMillis = System.currentTimeMillis(); break
        }

    }

    boolean isNew() { return status == NEW }

    boolean isSubmitted() { return status == SUBMITTED }

    boolean isRunning() { return status == RUNNING }

    boolean isCompleted()  { return status == COMPLETED  }

    boolean isActive() { status == SUBMITTED || status == RUNNING }

    protected StringBuilder toStringBuilder(StringBuilder builder) {
        builder << "id: ${task.id}; name: ${task.name}; status: $status; exit: ${task.exitStatus != Integer.MAX_VALUE ? task.exitStatus : '-'}; error: ${task.error ?: '-'}; workDir: ${task.workDir?.toUriString()}"
    }

    String toString() {
        def builder = toStringBuilder( new StringBuilder() )
        return "TaskHandler[${builder.toString()}]"
    }

    /**
     * The task status string. It extends the {@link TaskStatus} semantic adding specific status code string for
     * failed executions
     *
     * @return
     *      Can be either:
     *      - NEW: task has just been created and not yet submitted for execution
     *      - SUBMITTED: task has been submitted for execution
     *      - RUNNING: task is currently running
     *      - COMPLETED: task execution successfully completed
     *      - FAILED: task execution returned an error code
     *      - ABORTED: task execution was aborted by NF (likely because another task forced the workflow termination)
     */
    String getStatusString() {
        if( task.failed ) return 'FAILED'
        if( task.aborted ) return 'ABORTED'
        return this.status.toString()
    }

    TraceRecord safeTraceRecord() {
        try {
            return getTraceRecord()
        }
        catch (Exception e) {
            log.debug "Unable to get task trace record -- cause: ${e.message}", e
            return null
        }
    }

    /**
     * @return An {@link TraceRecord} instance holding task runtime information
     */
    @CompileDynamic
    TraceRecord getTraceRecord() {
        def record = new TraceRecord()
        record.task_id = task.id
        record.status = getStatusString()
        record.hash = task.hashLog
        record.name = task.name
        record.exit = task.exitStatus
        record.submit = this.submitTimeMillis
        record.start = this.startTimeMillis
        record.process = task.processor.getName()
        record.tag = task.config.tag
        record.module = task.config.module
        record.container = task.isContainerEnabled() ? task.getContainer() : null
        record.attempt = task.config.attempt

        record.script = task.getTraceScript()
        record.scratch = task.getScratch()
        record.workdir = task.getWorkDirStr()
        record.queue = task.config.queue
        record.cpus = task.config.getCpus()
        record.memory = task.config.getMemory()?.toBytes()
        record.disk = task.config.getDisk()?.toBytes()
        record.time = task.config.getTime()?.toMillis()
        record.env = task.getEnvironmentStr()
        record.executorName = task.processor.executor.getName()
        record.containerMeta = task.containerMeta()
        record.accelerator = task.config.getAccelerator()?.request
        record.accelerator_type = task.config.getAccelerator()?.type

        if( isCompleted() ) {
            record.error_action = task.errorAction?.toString()

            if( completeTimeMillis ) {
                // completion timestamp
                record.complete = completeTimeMillis
                // elapsed time since submit until completion
                if( submitTimeMillis )
                    record.duration = completeTimeMillis - submitTimeMillis
                // elapsed time since start of the job until completion
                // note: this may be override run time provided by the trace file (3rd line)
                if( startTimeMillis ) {
                    record.realtime = completeTimeMillis - startTimeMillis
                    log.trace "task stats: ${task.name}; start: ${startTimeMillis}; complete: ${completeTimeMillis}; realtime: ${completeTimeMillis - startTimeMillis} [${record.realtime}]; "
                }
            }

            // When Fusion trace is enabled, task metrics are collected by the Fusion client
            // and reported in .fusion/trace.json — skip parsing the .command.trace file
            // generated by the bash trace wrapper, as it won't exist
            final fusionTraceEnabled = task.processor.executor.isFusionEnabled() && Global.isFusionTraceEnabled()
            if( !fusionTraceEnabled ) {
                final file = task.workDir?.resolve(TaskRun.CMD_TRACE)
                try {
                    if(file) record.parseTraceFile(file)
                }
                catch( NoSuchFileException e ) {
                    log.trace "Unable to find trace file: $file"
                }
                catch( IOException e ) {
                    log.debug "[WARN] Cannot read trace file: $file -- Cause: ${e.message}"
                }
            }

            parseFusionTrace(record)
        }

        return record
    }

    protected void parseFusionTrace(TraceRecord record) {
        if( !task.processor.executor.isFusionEnabled() )
            return
        // Resolve the Fusion trace file path (.fusion/trace.json) in the task work directory.
        // This file is produced by the Fusion client and contains process, cgroup, and GPU metrics.
        final fusionTrace = task.workDir?.resolve(TaskRun.FUSION_TRACE)
        try {
            if( !fusionTrace )
                return
            // Parse the full Fusion trace JSON — sections: 'proc', 'cgroup', 'gpu'
            final json = TraceRecord.parseFusionTraceFile(fusionTrace)
            // GPU metrics are always extracted when available, regardless of the
            // NXF_FUSION_TRACE setting — this preserves backward compatibility
            final gpu = (Map<String,Object>) json.get('gpu')
            if( gpu )
                record.gpuMetrics = gpu
            // When NXF_FUSION_TRACE is disabled, stop here — only GPU metrics are collected.
            // When enabled, also map the 'proc' and 'cgroup' sections into TraceRecord fields,
            // replacing the metrics that would normally come from the bash trace wrapper
            if( !Global.isFusionTraceEnabled() )
                return
            applyFusionMetrics(record, json)
        }
        catch( NoSuchFileException e ) {
            log.trace "Unable to find Fusion trace file: $fusionTrace"
        }
        catch( Exception e ) {
            log.debug "[WARN] Cannot read Fusion trace file: $fusionTrace -- Cause: ${e.message}"
        }
    }

    protected void applyFusionMetrics(TraceRecord record, Map<String,Object> json) {
        final proc = (Map<String,Object>) json.get('proc')
        final cgroup = (Map<String,Object>) json.get('cgroup')

        if( proc ) {
            // CPU and timing metrics
            record.store.put('realtime', toLong(proc.get('realtime')))
            record.store.put('%cpu', toFloat(proc.get('pct_cpu')) / 10.0f as float)
            record.store.put('cpu_model', proc.get('cpu_name')?.toString())
            // I/O metrics
            record.store.put('rchar', toLong(proc.get('rchar')))
            record.store.put('wchar', toLong(proc.get('wchar')))
            record.store.put('syscr', toLong(proc.get('syscr')))
            record.store.put('syscw', toLong(proc.get('syscw')))
            record.store.put('read_bytes', toLong(proc.get('read_bytes')))
            record.store.put('write_bytes', toLong(proc.get('write_bytes')))
            // context switches
            record.store.put('vol_ctxt', toLong(proc.get('vol_ctxt')))
            record.store.put('inv_ctxt', toLong(proc.get('inv_ctxt')))
        }

        // Prefer cgroup memory metrics (more accurate for containerized tasks)
        if( cgroup ) {
            record.store.put('vmem', toLong(cgroup.get('memory_current')))
            record.store.put('rss', toLong(cgroup.get('memory_rss')))
            record.store.put('peak_vmem', toLong(cgroup.get('memory_peak')))
            record.store.put('peak_rss', toLong(cgroup.get('memory_peak_rss')))
            // Compute %mem as peak RSS against the cgroup memory limit — i.e. how close the
            // task got to its allocated memory. This overrides proc.pct_mem, which is RSS against
            // total host memory (ps %MEM semantics) and underestimates utilization for containerized
            // tasks whose limit is smaller than the host. Using peak (not current) RSS because Fusion
            // overwrites memory_rss post-exit, making the last sample unrepresentative.
            final memLimit = toLong(cgroup.get('memory_limit'))
            final memPeakRss = toLong(cgroup.get('memory_peak_rss'))
            if( memLimit > 0 ) {
                record.store.put('%mem', (memPeakRss / (float)memLimit * 100.0f) as float)
            }
        }
        else if( proc ) {
            // Fallback to proc memory metrics when cgroup is not available
            record.store.put('vmem', toLong(proc.get('vmem')))
            record.store.put('rss', toLong(proc.get('rss')))
            record.store.put('peak_vmem', toLong(proc.get('peak_vmem')))
            record.store.put('peak_rss', toLong(proc.get('peak_rss')))
            record.store.put('%mem', toFloat(proc.get('pct_mem')) / 10.0f as float)
        }
    }

    /**
     * Determine if a process can be forked i.e. can launch
     * a parallel task execution. This is only enforced when
     * when the `process.maxForks` directive is greater than zero.
     *
     * @return
     *      {@code true} if the number of forked process is less
     *      then of {@code process.maxForks} or {@code process.maxForks} is zero.
     *      {@code false} otherwise
     */
    boolean canForkProcess() {
        final max = task.processor.maxForks
        return !max ? true : task.processor.forksCount < max
    }

    /**
     * Determine if a task is ready for execution or it depends on resources
     * e.g. container that needs to be provisioned
     *
     * @return {@code true} when the task is ready for execution, {@code false} otherwise
     */
    boolean isReady() {
        task.isContainerReady()
    }

    /**
     * Increment the number of current forked processes
     */
    final void incProcessForks() {
        task.processor.forksCount?.increment()
    }

    /**
     * Decrement the number of current forked processes
     */
    final void decProcessForks() {
        task.processor.forksCount?.decrement()
    }

    /**
     * Check if the task submit could not be accomplished with the time specified via the
     * `maxWait` directive
     *
     * @return
     *      {@code true} if the task is in `submit` status after the amount of time specified
     *      via {@code maxAwait} directive has passed, otherwise {@code false} is returned.
     */
    boolean isSubmitTimeout() {
        final maxAwait = task.config.getMaxSubmitAwait()
        if( !maxAwait )
            return false
        final now = System.currentTimeMillis()
        if( isSubmitted() && now-submitTimeMillis>maxAwait.millis )
            return true
        return false
    }

    /**
     * Prepend the workflow Id to the job/task name. The workflow id is defined
     * by the environment variable {@code TOWER_WORKFLOW_ID}
     *
     * @param name
     *      The desired job name
     * @param env
     *      A map representing the variables in the host environment
     * @return
     *  The job having the prefix {@code tw-<ID>} when the variable {@code TOWER_WORKFLOW_ID}
     *  is defined in the host environment or just {@code name} otherwise
     */
    static String prependWorkflowPrefix(String name, Map<String,String> env) {
        final workflowId = env.get("TOWER_WORKFLOW_ID")
        return workflowId ? "tw-${workflowId}-${name}" : name
    }

}
