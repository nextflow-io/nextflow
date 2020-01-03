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

import groovy.text.GStringTemplateEngine
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import nextflow.script.WorkflowMetadata
/**
 * Render pipeline report processes execution.
 * Based on original TimelineObserver code by Paolo Di Tommaso
 *
 * @author Phil Ewels <phil.ewels@scilifelab.se>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ReportObserver implements TraceObserver {

    static final public String DEF_FILE_NAME = 'report.html'

    static final public int DEF_MAX_TASKS = 10_000

    /**
     * Holds the the start time for tasks started/submitted but not yet completed
     */
    final private Map<TaskId,TraceRecord> records = new LinkedHashMap<>()

    /**
     * Holds workflow session
     */
    private Session session

    /**
     * The path the HTML report file created
     */
    private Path reportFile

    /**
     * Max number of tasks allowed in the report, when they exceed this
     * number the tasks table is omitted
     */
    private int maxTasks = DEF_MAX_TASKS

    /**
     * Compute resources usage stats
     */
    private ResourcesAggregator aggregator

    /**
     * Creates a report observer
     *
     * @param file The file path where to store the resulting HTML report document
     */
    ReportObserver( Path file ) {
        this.reportFile = file
    }

    /**
     * Enables the collection of the task executions metrics in order to be reported in the HTML report
     *
     * @return {@code true}
     */
    @Override
    boolean enableMetrics() {
        return true
    }

    /**
     * @return The {@link WorkflowMetadata} object associated to this execution
     */
    protected WorkflowMetadata getWorkflowMetadata() {
        session.getWorkflowMetadata()
    }

    /**
     * @return The map of collected {@link TraceRecord}s
     */
    protected Map<TaskId,TraceRecord> getRecords() {
        records
    }

    /**
     * Set the number max allowed tasks. If this number is exceed the the tasks
     * json in not included in the final report
     *
     * @param value The number of max task record allowed to be included in the HTML report
     * @return The {@link ReportObserver} itself
     */
    ReportObserver setMaxTasks( int value ) {
        this.maxTasks = value
        return this
    }

    /**
     * Create the trace file, in file already existing with the same name it is
     * "rolled" to a new file
     */
    @Override
    void onFlowCreate(Session session) {
        this.session = session
        this.aggregator = new ResourcesAggregator(session)
    }

    /**
     * Save the pending processes and close the trace file
     */
    @Override
    void onFlowComplete() {
        log.debug "Flow completing -- rendering html report"
        try {
            renderHtml()
        }
        catch (Exception e) {
            log.warn "Failed to render execution report -- see the log file for details", e
        }
    }

    /**
     * This method is invoked when a process is created
     *
     * @param process A {@link TaskProcessor} object representing the process created
     */
    @Override
    void onProcessCreate(TaskProcessor process) { }


    /**
     * This method is invoked before a process run is going to be submitted
     *
     * @param handler A {@link TaskHandler} object representing the task submitted
     */
    @Override
    void onProcessSubmit(TaskHandler handler, TraceRecord trace) {
        log.trace "Trace report - submit process > $handler"
        synchronized (records) {
            records[ trace.taskId ] = trace
        }
    }

    /**
     * This method is invoked when a process run is going to start
     *
     * @param handler  A {@link TaskHandler} object representing the task started
     */
    @Override
    void onProcessStart(TaskHandler handler, TraceRecord trace) {
        log.trace "Trace report - start process > $handler"
        synchronized (records) {
            records[ trace.taskId ] = trace
        }
    }

    /**
     * This method is invoked when a process run completes
     *
     * @param handler A {@link TaskHandler} object representing the task completed
     */
    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        log.trace "Trace report - complete process > $handler"
        if( !trace ) {
            log.debug "WARN: Unable to find trace record for task id=${handler.task?.id}"
            return
        }

        synchronized (records) {
            records[ trace.taskId ] = trace
            aggregate(trace)
        }
    }

    /**
     * This method is invoked when a process run cache is hit
     *
     * @param handler A {@link TaskHandler} object representing the task cached
     */
    @Override
    void onProcessCached(TaskHandler handler, TraceRecord trace) {
        log.trace "Trace report - cached process > $handler"

        // event was triggered by a stored task, ignore it
        if( trace == null ) {
            return
        }

        // remove the record from the current records
        synchronized (records) {
            records[ trace.taskId ] = trace
            aggregate(trace)
        }
    }

    /**
     * Aggregates task record for each process in order to render the
     * final execution stats
     *
     * @param record A {@link TraceRecord} object representing a task executed
     */
    protected void aggregate(TraceRecord record) {
        aggregator.aggregate(record)
    }

    /**
     * @return The tasks json payload
     */
    protected String renderTasksJson() {
        final r = getRecords()
        r.size()<=maxTasks ? renderJsonData(r.values()) : 'null'
    }

    protected String renderSummaryJson() {
        final result = aggregator.renderSummaryJson()
        log.debug "Execution report summary data:\n  ${result}"
        return result
    }

    protected String renderPayloadJson() {
        "{ \"trace\":${renderTasksJson()}, \"summary\":${renderSummaryJson()} }"
    }

    /**
     * Render the report HTML document
     */
    protected void renderHtml() {

        // render HTML report template
        final tpl_fields = [
            workflow : getWorkflowMetadata(),
            payload : renderPayloadJson(),
            assets_css : [
                readTemplate('assets/bootstrap.min.css'),
                readTemplate('assets/datatables.min.css')
            ],
            assets_js : [
                readTemplate('assets/jquery-3.2.1.min.js'),
                readTemplate('assets/popper.min.js'),
                readTemplate('assets/bootstrap.min.js'),
                readTemplate('assets/datatables.min.js'),
                readTemplate('assets/moment.min.js'),
                readTemplate('assets/plotly.min.js'),
                readTemplate('assets/ReportTemplate.js')
            ]
        ]
        final tpl = readTemplate('ReportTemplate.html')
        def engine = new GStringTemplateEngine()
        def html_template = engine.createTemplate(tpl)
        def html_output = html_template.make(tpl_fields).toString()

        // make sure the parent path exists
        def parent = reportFile.getParent()
        if( parent )
            Files.createDirectories(parent)

        // roll the any trace files that may exist
        reportFile.rollFile()

        def writer = Files.newBufferedWriter(reportFile, Charset.defaultCharset())
        writer.withWriter { w -> w << html_output }
        writer.close()
    }

    /**
     * Render the executed tasks json payload
     *
     * @param data A collection of {@link TraceRecord}s representing the tasks executed
     * @return The rendered json payload
     */
    protected String renderJsonData(Collection<TraceRecord> data) {
        def List<String> formats = null
        def List<String> fields = null
        def result = new StringBuilder()
        result << '[\n'
        int i=0
        for( TraceRecord record : data ) {
            if( i++ ) result << ','
            if( !formats ) formats = TraceRecord.FIELDS.values().collect { it!='str' ? 'num' : 'str' }
            if( !fields ) fields = TraceRecord.FIELDS.keySet() as List
            record.renderJson(result,fields,formats)
        }
        result << ']'
        return result.toString()
    }

    /**
     * Read the document HTML template from the application classpath
     *
     * @param path A resource path location
     * @return The loaded template as a string
     */
    private String readTemplate( String path ) {
        StringWriter writer = new StringWriter();
        def res =  this.class.getResourceAsStream( path )
        int ch
        while( (ch=res.read()) != -1 ) {
            writer.append(ch as char);
        }
        writer.toString();
    }

}
