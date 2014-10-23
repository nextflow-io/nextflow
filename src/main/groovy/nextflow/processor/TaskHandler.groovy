/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
import java.util.concurrent.CountDownLatch

import groovy.util.logging.Slf4j
import nextflow.trace.TraceRecord
/**
 * Actions to handle the underlying job running the user task.
 *
 * <p>
 * Note this model the job in the execution facility (i.e. grid, cloud, etc)
 * NOT the *logical* user task
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
public abstract class TaskHandler {

    enum Status { NEW, SUBMITTED, RUNNING, COMPLETED }

    protected TaskHandler(TaskRun task, TaskConfig taskConfig) {
        this.task = task
        this.taskConfig = taskConfig
    }


    /** Only for testing purpose */
    protected TaskHandler() {

    }

    /**
     * The task managed by this handler
     */
    TaskRun task

    /**
     * The configuration object defined by this task
     */
    TaskConfig taskConfig

    /**
     * The task managed by this handler
     */
    TaskRun getTask() { task }

    /**
     * Task current status
     */
    volatile Status status = Status.NEW

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
     * @param status The sask status as defined by {@link Status}
     */
    def void setStatus( Status status ) {

        // skip if the status is the same aam
        if ( this.status == status || status == null )
            return

        // change the status
        this.status = status
        switch( status ) {
            case Status.SUBMITTED: submitTimeMillis = System.currentTimeMillis(); break
            case Status.RUNNING: startTimeMillis = System.currentTimeMillis(); break
            case Status.COMPLETED: completeTimeMillis = System.currentTimeMillis(); break
        }

    }

    boolean isNew() { return status == Status.NEW }

    boolean isSubmitted() { return status == Status.SUBMITTED }

    boolean isRunning() { return status == Status.RUNNING }

    boolean isCompleted()  { return status == Status.COMPLETED  }

    protected StringBuilder toStringBuilder(StringBuilder builder) {
        builder << "id: ${task.id}; name: ${task.name}; status: $status; exit: ${task.exitStatus != Integer.MAX_VALUE ? task.exitStatus : '-'}; workDir: ${task.workDir}"
    }

    String toString() {
        def builder = toStringBuilder( new StringBuilder() )
        return "TaskHandler[${builder.toString()}]"
    }

    /**
     * @return An {@link TraceRecord} instance holding task runtime information
     */
    TraceRecord getTraceRecord() {
        def record = new TraceRecord()
        record.task_id = task.id
        record.status = task.failed ? 'FAILED' : this.status
        record.hash = task.hashLog
        record.name = task.name
        record.exit_status = task.exitStatus
        record.submit = this.submitTimeMillis
        record.start = this.startTimeMillis

        if( isCompleted() ) {
            if( completeTimeMillis ) {
                // completion timestamp
                record.complete = completeTimeMillis
                // elapsed time since submit until completion
                if( submitTimeMillis )
                    record.wall_time = completeTimeMillis - submitTimeMillis
                // elapsed time since start of the job until completion
                // note: this may be override run time provided by the trace file (3rd line)
                if( startTimeMillis )
                    record.run_time = completeTimeMillis - startTimeMillis
            }

            def file = task.workDir?.resolve(TaskRun.CMD_TRACE)
            if( file?.exists() ) {
                record.parseTraceFile(file.text)
            }
        }

        return record
    }


}
