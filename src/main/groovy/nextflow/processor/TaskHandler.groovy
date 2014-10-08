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

    TraceRecord getTraceRecord() {
        def record = new TraceRecord()
        record.task_id = task.id
        record.status = this.status
        record.hash = task.hashLog
        record.name = task.name
        record.exit_status = task.exitStatus
        record.submit = this.submitTimeMillis
        record.start = this.startTimeMillis
        record.complete = this.completeTimeMillis

        if( isCompleted() ) {
            def file = task.workDir?.resolve(TaskRun.CMD_TRACE)
            log.debug "Command trace file exists: ${file.exists()} -- $file"

            if( file?.exists() ) {
                record.putAll( parseTraceFile(file.text) )
            }
        }

        return record
    }

    private int traceParseError

    private static String[] VALID = 'rchar,wchar,syscr,syscw,read_bytes,write_bytes,Threads'.split(',')

    /**
     * parse the trace file
     */

    Map<String,Object> parseTraceFile( String text ) {

        String[] header = null
        final result = [:]
        final lines = text.readLines()
        for( int line=0; line<lines.size(); line++ ) {
            String row = lines[line]

            /*
             * 1st line -- parse the header
             */
            if( line == 0 ) {
                header = row.trim().split(/\s+/)
            }
            /*
             * 2nd line -- parse values produced by 'ps'
             */
            else if( line == 1 ) {
                String[] values = row.trim().split(/\s+/)
                if( header.length == 5 && values.length == 5 ) {
                    for( int i=0; i<values.length; i++ ) {
                        def name = header[i]
                        if( name == 'rss' || name == 'vmem')
                            result[name] = values[i].toLong() * 1024
                        else {
                            result[name] = values[i]
                        }
                    }
                }
                else {
                    if( traceParseError++ == 0 ) {
                        log.debug "Invalid process file format: \n${text.indent('  ')}\n"
                    }
                    continue
                }
            }
            /*
             * 3rd line -- max of the above values
             */
            else if( line == 2 ) {
                String[] values = row.trim().split(/\s+/)
                for( int i=0; i<values.length; i++ ) {
                    def name = header[i]
                    if( name == 'rss' || name == 'vmem')
                        result["max_$name"] = values[i].toLong() * 1024
                    else {
                        result["max_$name"] = values[i]
                    }
                }

            }
            /*
             * all remaining are in the form <key: value>
             */
            else {
                def items = row.trim().split(/:\s+/)
                if( !items ) continue

                if( items.size() == 2 && items[0] in TaskHandler.VALID )  {
                    result[items[0].toLowerCase()] = items[1]
                }
                else if ( items[0].startsWith('Vm')  ) {
                    def val = items[1]
                    if( val.endsWith(' kB') ) {
                        val = items[1][0..-3].toLong() * 1024
                    }

                    result[items[0].toLowerCase()] = val
                }

            }
        }

        return result
    }


}
