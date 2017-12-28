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

    static final String DEF_FILE_NAME = 'report.html'


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

    ReportObserver( Path file ) {
        this.reportFile = file
    }

    /**
     * @return The {@link WorkflowMetadata} object associated to this execution
     */
    protected WorkflowMetadata getWorkflowMetadata() {
        session.binding.getVariable('workflow') as WorkflowMetadata
    }

    protected Map<TaskId,TraceRecord> getRecords() { records }

    /**
     * Create the trace file, in file already existing with the same name it is
     * "rolled" to a new file
     */
    @Override
    void onFlowStart(Session session) {
        this.session = session
    }

    /**
     * Save the pending processes and close the trace file
     */
    @Override
    void onFlowComplete() {
        log.debug "Flow completing -- rendering html report"
        renderHtml()
    }


    @Override
    void onProcessCreate(TaskProcessor process) { }


    /**
     * This method is invoked before a process run is going to be submitted
     * @param handler
     */
    @Override
    void onProcessSubmit(TaskHandler handler) {
        log.trace "Trace report - submit process > $handler"
        final trace = handler.getTraceRecord()
        synchronized (records) {
            records[ trace.taskId ] = trace
        }
    }

    /**
     * This method is invoked when a process run is going to start
     * @param handler
     */
    @Override
    void onProcessStart(TaskHandler handler) {
        log.trace "Trace report - start process > $handler"
        def trace = handler.getTraceRecord()
        synchronized (records) {
            records[ trace.taskId ] = trace
        }
    }

    /**
     * This method is invoked when a process run completes
     * @param handler
     */
    @Override
    void onProcessComplete(TaskHandler handler) {
        log.trace "Trace report - complete process > $handler"
        final record = handler.getTraceRecord()
        if( !record ) {
            log.debug "WARN: Unable to find trace record for task id=${handler.task?.id}"
            return
        }

        // remove the record from the current records
        synchronized (records) {
            records[ record.taskId ] = record
        }
    }

    @Override
    void onProcessCached(TaskHandler handler) {
        log.trace "Trace report - cached process > $handler"
        final record = handler.getTraceRecord()

        // remove the record from the current records
        synchronized (records) {
            records[ record.taskId ] = record
        }

    }

    protected void renderHtml() {

        // render HTML report template
        final tpl_fields = [
            workflow : getWorkflowMetadata(),
            payload : renderJsonData(records.values()),
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

    protected String renderJsonData(Collection<TraceRecord> data) {
        def List<String> formats = null
        def List<String> fields = null
        def result = new StringBuilder()
        result << '{ "trace": [\n'
        int i=0
        for( TraceRecord record : data ) {
            if( i++ ) result << ','
            if( !formats ) formats = TraceRecord.FIELDS.values().collect { it!='str' ? 'num' : 'str' }
            if( !fields ) fields = TraceRecord.FIELDS.keySet() as List
            record.renderJson(result,fields,formats)
        }
        result << ']}'
        return result.toString()
    }

    protected String readTemplate( String path ) {
        StringWriter writer = new StringWriter();
        def res =  this.class.getResourceAsStream( path )
        int ch
        while( (ch=res.read()) != -1 ) {
            writer.append(ch as char);
        }
        writer.toString();
    }

    @Override
    boolean enableMetrics() {
        return true
    }
}
