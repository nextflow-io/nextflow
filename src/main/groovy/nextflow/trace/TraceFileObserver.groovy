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
    final static public List<String> HEAD = [
            'task-id',
            'hash',
            'native-id',
            'name',
            'status',
            'exit-status',
            'submit',
            'start',
            'complete',
            'wall-time',
            'run-time',
            'mem',
            'cpu',
            'io'
    ]

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
    private Map<Object,TraceRecord> current = new ConcurrentHashMap<>()

    private Agent<PrintWriter> writer


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
        current.values().each { record -> traceFile.println(render(record)) }
        traceFile.flush()
        traceFile.close()
    }


    /**
     * This method is invoked before a process run is going to be submitted
     * @param handler
     */
    void onProcessSubmit(TaskHandler handler) {
        def trace = handler.getTraceRecord()
        current[ trace.taskId ] = trace
    }

    /**
     * This method is invoked when a process run is going to start
     * @param handler
     */
    void onProcessStart(TaskHandler handler) {
        def trace = handler.getTraceRecord()
        current[ trace.taskId ] = trace
    }

    /**
     * This method is invoked when a process run completes
     * @param handler
     */
    void onProcessComplete(TaskHandler handler) {
        final taskId = handler.task.id
        final record = handler.getTraceRecord()
        if( !record ) {
            log.debug "Profile warn: Unable to find record for task_run with id: ${taskId}"
            return
        }

        // remove the record from the current records
        current.remove(taskId)

        // save to the file
        writer.send { PrintWriter it -> it.println(render(record)) }
    }

    /**
     * Render a {@link TraceRecord} object to a string
     *
     * @param trace
     * @return
     */
    String render(TraceRecord trace) {
        assert trace
        def all = new ArrayList(14)
        trace.with{
            def wallTime = complete && submit ? complete-submit : 0
            def runTime = complete && start ? complete-start : 0

            all << taskId
            all << str(hash)
            all << str(nativeId)
            all << str(name)
            all << str(status)
            all << str(exit)
            all << date(submit)
            all << date(start)
            all << date(complete)
            all << wallTime
            all << runTime
            all << str(mem)
            all << str(cpu)
            all << str(io)
        }

        all.join(delim)
    }

    static String str(item) {
        return item ? item.toString() : '-'
    }

    static String date(long item) {
        if( !item ) return '-'
        FMT.format(new Date(item))
    }

}
