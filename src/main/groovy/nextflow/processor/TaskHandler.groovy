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

package nextflow.processor

import static nextflow.processor.TaskStatus.*

import java.nio.file.NoSuchFileException
import java.util.concurrent.CountDownLatch

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
public abstract class TaskHandler {

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

    CountDownLatch latch

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
     * Force the submitted job to quit
     */
    abstract void kill()

    /**
     * Submit the task for execution.
     *
     * Note: the underlying execution platform may schedule it in its own queue
     */
    abstract void submit()

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

    /**
     * @return An {@link TraceRecord} instance holding task runtime information
     */
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
        record.container = task.container
        record.attempt = task.config.attempt

        record.script = task.getScript()
        record.scratch = task.getScratch()
        record.workdir = task.getWorkDirStr()
        record.queue = task.config.queue
        record.cpus = task.config.getCpus()
        record.memory = task.config.getMemory()?.toBytes()
        record.disk = task.config.getDisk()?.toBytes()
        record.time = task.config.getTime()?.toMillis()
        record.env = task.getEnvironmentStr()

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

            def file = task.workDir?.resolve(TaskRun.CMD_TRACE)
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


}
