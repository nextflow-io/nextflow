/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
            'peak_rss',
            'peak_vmem',
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
    void onProcessSubmit(TaskHandler handler, TraceRecord trace) {
        current[ trace.taskId ] = trace
    }

    /**
     * This method is invoked when a process run is going to start
     * @param handler
     */
    @Override
    void onProcessStart(TaskHandler handler, TraceRecord trace) {
        current[ trace.taskId ] = trace
    }

    /**
     * This method is invoked when a process run completes
     * @param handler
     */
    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        final taskId = handler.task.id
        if( !trace ) {
            log.debug "Profile warn: Unable to find record for task_run with id: ${taskId}"
            return
        }

        // remove the record from the current records
        current.remove(taskId)

        // save to the file
        writer.send { PrintWriter it -> it.println(render(trace)); it.flush() }
    }

    @Override
    void onProcessCached(TaskHandler handler, TraceRecord trace) {
        // save to the file
        writer.send { PrintWriter it -> it.println(render( trace )); it.flush() }
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
