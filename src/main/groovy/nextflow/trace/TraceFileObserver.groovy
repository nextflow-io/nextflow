/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
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
import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardOpenOption
import java.util.concurrent.ConcurrentHashMap

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.agent.Agent
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
/**
 * Create a CSV file containing the processes execution information
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TraceFileObserver implements TraceObserver {

    static final String DEF_FILE_NAME = 'trace.txt'

    /**
     * The list of fields included in the trace report
     */
    List<String> fields = [
            'task_id',
            'hash',
            'native_id',
            'name',
            'status',
            'exit',
            'submit',
            'duration',
            'realtime',
            '%cpu',
            'rss',
            'vmem',
            'rchar',
            'wchar'
    ]

    List<String> formats


    /**
     * The delimiter character used to separate column in the CSV file
     */
    String separator = '\t'

    /**
     * The path where the file is created. It is set by the object constructor
     */
    private Path tracePath

    /**
     * The actual file object
     */
    private PrintWriter traceFile

    /**
     * Holds the the start time for tasks started/submitted but not yet completed
     */
    @PackageScope Map<TaskId,TraceRecord> current = new ConcurrentHashMap<>()

    private Agent<PrintWriter> writer

    private boolean useRawNumber

    void setFields( List<String> entries ) {

        def names = TraceRecord.FIELDS.keySet()
        def result = new ArrayList<String>(entries.size())
        for( def item : entries ) {
            def thisName = item.trim()

            if( thisName ) {
                if( thisName in names )
                    result << thisName
                else {
                    String message = "Not a valid trace field name: '$thisName'"
                    def alternatives = names.bestMatches(thisName)
                    if( alternatives )
                        message += " -- Possible solutions: ${alternatives.join(', ')}"
                    throw new IllegalArgumentException(message)
                }
            }

        }

        this.fields = result
    }

    TraceFileObserver setFieldsAndFormats( value ) {
        List<String> entries
        if( value instanceof String ) {
            entries = value.tokenize(', ')
        }
        else if( value instanceof List ) {
            entries = (List)value
        }
        else {
            throw new IllegalArgumentException("Not a valid trace fields value: $value")
        }

        List<String> fields = new ArrayList<>(entries.size())
        List<String> formats = new ArrayList<>(entries.size())

        for( String x : entries ) {
            String name
            String fmt
            int p = x.indexOf(':')
            if( p == -1 ) {
                name = x
                fmt = TraceRecord.FIELDS.get(name)      // get the default type
            }
            else {
                name = x.substring(0,p)
                fmt = x.substring(p+1)
            }

            if( !fmt )
                throw new IllegalArgumentException("Unknown trace field name: `$name`")

            if( useRawNumber && fmt in TraceRecord.NON_PRIMITIVE_TYPES ) {
                fmt = 'num'
            }

            fields << name.trim()
            formats << fmt.trim()
        }

        setFields(fields)
        setFormats(formats)

        return this
    }

    TraceFileObserver useRawNumbers( boolean value ) {
        this.useRawNumber = value

        List<String> local = []
        for( String name : fields ) {
            def type = TraceRecord.FIELDS.get(name)
            if( useRawNumber && type in TraceRecord.NON_PRIMITIVE_TYPES ) {
                type = 'num'
            }
            local << type
        }
        this.formats = local
        return this
    }


    /**
     * Create the trace observer
     *
     * @param traceFile A path to the file where save the tracing data
     */
    TraceFileObserver( Path traceFile ) {
        this.tracePath = traceFile
    }

    /** ONLY FOR TESTING PURPOSE */
    protected TraceFileObserver( ) {}

    /**
     * Create the trace file, in file already existing with the same name it is
     * "rolled" to a new file
     */
    @Override
    void onFlowStart(Session session) {
        log.debug "Flow starting -- trace file: $tracePath"

        // make sure parent path exists
        def parent = tracePath.getParent()
        if( parent )
            Files.createDirectories(parent)

        // roll the any trace files that may exist
        tracePath.rollFile()

        // create a new trace file
        traceFile = new PrintWriter(Files.newBufferedWriter(tracePath, Charset.defaultCharset(), StandardOpenOption.APPEND, StandardOpenOption.CREATE))

        // launch the agent
        writer = new Agent<PrintWriter>(traceFile)
        writer.send { traceFile.println(fields.join(separator)); traceFile.flush() }
    }

    /**
     * Save the pending processes and close the trace file
     */
    @Override
    void onFlowComplete() {
        log.debug "Flow completing -- flushing trace file"
        // wait for termination and flush the agent content
        writer.await()

        // write the remaining records
        current.values().each { record -> traceFile.println(render(record)) }
        traceFile.flush()
        traceFile.close()
    }


    @Override
    void onProcessCreate(TaskProcessor process) {

    }


    /**
     * This method is invoked before a process run is going to be submitted
     * @param handler
     */
    @Override
    void onProcessSubmit(TaskHandler handler) {
        def trace = handler.getTraceRecord()
        current[ trace.taskId ] = trace
    }

    /**
     * This method is invoked when a process run is going to start
     * @param handler
     */
    @Override
    void onProcessStart(TaskHandler handler) {
        def trace = handler.getTraceRecord()
        current[ trace.taskId ] = trace
    }

    /**
     * This method is invoked when a process run completes
     * @param handler
     */
    @Override
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
        writer.send { PrintWriter it -> it.println(render(record)); it.flush() }
    }

    @Override
    void onProcessCached(TaskHandler handler) {
        // save to the file
        writer.send { PrintWriter it -> it.println(render( handler.getTraceRecord() )); it.flush() }
    }

    /**
     * Render a {@link TraceRecord} object to a string
     *
     * @param trace
     * @return
     */
    String render(TraceRecord trace) {
        assert trace
        trace.renderText(fields, formats, separator)
    }

    @Override
    boolean enableMetrics() {
        return true
    }
}
