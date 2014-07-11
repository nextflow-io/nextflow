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

package nextflow.trace
import java.nio.file.Path
import java.text.SimpleDateFormat
import java.util.concurrent.ConcurrentHashMap

import groovy.util.logging.Slf4j
import groovyx.gpars.agent.Agent
import nextflow.Session
import nextflow.extension.FilesExtensions
import nextflow.processor.TaskHandler
/**
 * Create a CSV file containing the processes execution information
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class TraceFileObserver implements TraceObserver {

    /**
     * The used date format to convert {@link Date} to strings
     */
    final static private SimpleDateFormat FMT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.SSS")

    /**
     * The CSV header
     */
    final static private List<String> HEAD = ['id','process name','submit','start','complete','wall-time','run-time']

    /**
     * The path where the file is created. It is set by the object constructor
     */
    private Path tracePath

    /**
     * The delimiter character used to separate column in the CSV file
     */
    private String delim = ','

    /**
     * The actual file object
     */
    private PrintWriter traceFile

    /**
     * Holds the the start time for tasks started/submitted but not yet completed
     */
    private Map<Object,ProfileRecord> current = new ConcurrentHashMap<>()

    private Agent<PrintWriter> writer

    /*
     * This object represent holds the information of a single process run,
     * its content is saved to a trace file line
     */
    class ProfileRecord {
        def id
        String name
        long submit
        long start
        long complete

        String toString() {
            def wallTime = complete && submit ? complete-submit : 0
            def runTime = complete && start ? complete-start : 0

            def builder = new StringBuilder()
            builder << id << delim
            builder << name << delim
            builder << fmt(submit) << delim
            builder << fmt(start) << delim
            builder << fmt(complete) << delim
            builder << wallTime << delim
            builder << runTime

            builder.toString()
        }

        static String fmt(long item) {
            if( !item ) return '-'
            FMT.format(new Date(item))
        }
    }


    /**
     * Create the trace observer
     *
     * @param traceFile A path to the file where save the tracing data
     */
    TraceFileObserver( Path traceFile ) {
        this.tracePath = traceFile
    }

    /**
     * Create the trace file, in file already existing with the same name it is
     * "rolled" to a new file
     */
    void onFlowStart(Session session) {
        log.debug "Flow starting -- trace file: $tracePath"
        // roll the any trace files that may exist
        FilesExtensions.rollFile(tracePath)

        // create a new trace file
        traceFile = new PrintWriter(new BufferedWriter( new FileWriter(tracePath.toFile())))
        traceFile.println( HEAD.join(delim) )

        // launch the agent
        writer = new Agent<PrintWriter>(traceFile)
    }

    /**
     * Save the pending processes and close the trace file
     */
    void onFlowComplete() {
        log.debug "Flow completing -- flushing trace file"
        // wait for termination and flush the agent content
        writer.await()

        // write the remaining records
        current.each { record -> traceFile.println(record) }
        traceFile.flush()
        traceFile.close()
    }


    /**
     * This method is invoked before a process run is going to be submitted
     * @param handler
     */
    void onProcessSubmit(TaskHandler handler) {
        def taskId = handler.task.id
        def taskName = handler.task.name
        def record = new ProfileRecord()
        record.with {
            id = taskId
            name = taskName
            submit = System.currentTimeMillis()
        }

        current[ taskId ] = record
    }

    /**
     * This method is invoked when a process run is going to start
     * @param handler
     */
    void onProcessStart(TaskHandler handler) {
        def record = current.get(handler.task.id)
        if( record ) {
            record.start = System.currentTimeMillis()
        }
        else {
            log.debug "Profile warn: Unable to find record for task_run with id: ${handler.task.id}"
        }
    }

    /**
     * This method is invoked when a process run completes
     * @param handler
     */
    void onProcessComplete(TaskHandler handler) {
        final taskId = handler.task.id
        final record = current.get(taskId)
        if( !record ) {
            log.debug "Profile warn: Unable to find record for task_run with id: ${taskId}"
            return
        }

        // update the record and remove it from the current records
        record.complete = System.currentTimeMillis()
        current.remove(taskId)
        // save to the file
        writer.send { PrintWriter it -> it.println(record) }
    }


}
