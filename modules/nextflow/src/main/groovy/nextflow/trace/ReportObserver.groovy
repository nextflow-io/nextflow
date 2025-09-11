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

package nextflow.trace

import java.nio.file.Files
import java.nio.file.Path

import groovy.text.GStringTemplateEngine
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.processor.TaskId
import nextflow.script.WorkflowMetadata
import nextflow.trace.config.ReportConfig
import nextflow.trace.event.TaskEvent
import nextflow.util.SysHelper
import nextflow.util.TestOnly
/**
 * Render pipeline report processes execution.
 * Based on original TimelineObserver code by Paolo Di Tommaso
 *
 * @author Phil Ewels <phil.ewels@scilifelab.se>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ReportObserver implements TraceObserverV2 {

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
    private int maxTasks

    /**
     * Compute resources usage stats
     */
    private ResourcesAggregator aggregator

    /**
     * Overwrite existing trace file (required in some cases, as rolling filename has been deprecated)
     */
    private boolean overwrite

    ReportObserver(ReportConfig config) {
        this.reportFile = FileHelper.asPath(config.file)
        this.maxTasks = config.maxTasks
        this.overwrite = config.overwrite
    }

    @TestOnly
    protected ReportObserver() {}

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
     * Create the trace file, in file already existing with the same name it is
     * "rolled" to a new file
     */
    @Override
    void onFlowCreate(Session session) {
        this.session = session
        this.aggregator = new ResourcesAggregator(session)
        // check if the process exists
        if( Files.exists(reportFile) && !overwrite )
            throw new AbortOperationException("Report file already exists: ${reportFile.toUriString()} -- enable the 'report.overwrite' option in your config file to overwrite existing files")
    }

    /**
     * Save the pending processes and close the trace file
     */
    @Override
    void onFlowComplete() {
        log.debug "Workflow completed -- rendering execution report"
        try {
            renderHtml()
        }
        catch (Exception e) {
            log.warn "Failed to render execution report -- see the log file for details", e
        }
    }

    @Override
    void onTaskSubmit(TaskEvent event) {
        log.trace "Trace report - submit process > $event.handler"
        synchronized (records) {
            records[ event.trace.taskId ] = event.trace
        }
    }

    @Override
    void onTaskStart(TaskEvent event) {
        log.trace "Trace report - start process > $event.handler"
        synchronized (records) {
            records[ event.trace.taskId ] = event.trace
        }
    }

    @Override
    void onTaskComplete(TaskEvent event) {
        log.trace "Trace report - complete process > $event.handler"
        if( !event.trace ) {
            log.debug "WARN: Unable to find trace record for task id=${event.handler.task?.id}"
            return
        }

        synchronized (records) {
            records[ event.trace.taskId ] = event.trace
            aggregate(event.trace)
        }
    }

    @Override
    void onTaskCached(TaskEvent event) {
        log.trace "Trace report - cached process > $event.handler"

        // event was triggered by a stored task, ignore it
        if( !event.trace ) {
            return
        }

        // remove the record from the current records
        synchronized (records) {
            records[ event.trace.taskId ] = event.trace
            aggregate(event.trace)
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
        aggregator.renderSummaryJson()
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
            // Add SysHelper to template binding so it can be used in templates
            SysHelper : SysHelper,
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
        final engine = new GStringTemplateEngine()
        final html_template = engine.createTemplate(tpl)
        final html_output = html_template.make(tpl_fields).toString()

        // make sure the parent path exists
        final parent = reportFile.getParent()
        if( parent )
            Files.createDirectories(parent)

        final writer = TraceHelper.newFileWriter(reportFile, overwrite, 'Report')
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
        List<String> formats = null
        List<String> fields = null
        final result = new StringBuilder()
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
        final writer = new StringWriter()
        final res = this.class.getResourceAsStream( path )
        int ch
        while( (ch=res.read()) != -1 ) {
            writer.append(ch as char)
        }
        writer.toString()
    }

}
