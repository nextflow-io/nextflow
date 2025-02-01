/*
 * Copyright 2013-2024, Seqera Labs
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

import java.nio.file.NoSuchFileException
import java.util.concurrent.atomic.AtomicBoolean

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.trace.TraceRecord
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

    /** Only for testing purpose */
    protected TaskHandler() { }

    /**
     * The task managed by this handler
     */
    TaskRun task

    /**
     * The task managed by this handler
     */
    TaskRun getTask() { task }

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
        record.container = task.getContainer()
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

            final file = task.workDir?.resolve(TaskRun.CMD_TRACE)
            try {
                if(file) record.parseTraceFile(file)
            }
            catch( NoSuchFileException e ) {
                // ignore it
            }
            catch( IOException e ) {
                log.debug "[WARN] Cannot read trace file: $file -- Cause: ${e.message}"
            }
        }

        return record
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
     * e.g. container that needs to be provisionied
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
