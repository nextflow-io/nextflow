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
import java.util.concurrent.ConcurrentHashMap

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import nextflow.util.Duration
import org.apache.commons.lang.StringEscapeUtils
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

    static final String DEF_FILE_NAME = 'nf-report.html'

    /**
     * Holds the the start time for tasks started/submitted but not yet completed
     */
    final private Map<TaskId,TraceRecord> records = new LinkedHashMap<>()

    /**
     * Holds workflow session
     */
    private Class nfsession

    // BREAKS: unable to resolve class WorkflowMetadata
    private WorkflowMetadata getWorkflowMetadata() {
        nfsession.binding.getVariable('workflow')
    }

    private Path reportFile

    private long beginMillis

    private long startMillis

    private long endMillis

    ReportObserver( Path file ) {
        this.reportFile = file
    }

    /**
     * Create the trace file, in file already existing with the same name it is
     * "rolled" to a new file
     */
    @Override
    void onFlowStart(Session session) {
        beginMillis = startMillis = System.currentTimeMillis()
        nfsession = session
    }

    /**
     * Save the pending processes and close the trace file
     */
    @Override
    void onFlowComplete() {
        log.debug "Flow completing -- rendering html report"
        endMillis = System.currentTimeMillis()
        // println session_binding
        // session_binding.each{
        //   println it.key
        //   println it.value
        //   println "-----"
        // }
        // workflow = session_binding.binding.getVariable('workflow')
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
        final taskId = handler.task.id
        final record = handler.getTraceRecord()
        if( !record ) {
            log.debug "Profile warn: Unable to find record for task_run with id: ${taskId}"
            return
        }

        // remove the record from the current records
        synchronized (records) {
            records[ record.taskId ] = record
        }
    }

    @Override
    void onProcessCached(TaskHandler handler) {

        final record = handler.getTraceRecord()

        // remove the record from the current records
        synchronized (records) {
            records[ record.taskId ] = record
        }

        beginMillis = Math.min( beginMillis, record.get('submit') as long )
    }

    protected void renderHtml() {

        // workflow metadata fields
        def metadata = getWorkflowMetadata()

        // make records safe for JS
        def records_safe = [:]
        records.values().eachWithIndex { TraceRecord it, index ->
            records_safe[index] = makesafe(it)
        }

        // render HTML report template
        final tpl_fields = [
            workflow : metadata,
            records : records_safe,
            assets_css : [
                readTemplate('assets/bootstrap.min.css'),
                readTemplate('assets/datatables.min.css')
            ],
            assets_js : [
                readTemplate('assets/jquery-3.2.1.min.js'),
                readTemplate('assets/popper.min.js'),
                readTemplate('assets/bootstrap.min.js'),
                readTemplate('assets/datatables.min.js'),
                readTemplate('assets/plotly.min.js'),
                readTemplate('assets/ReportTemplate.js')
            ]
        ]
        final tpl = readTemplate('ReportTemplate.html')
        def engine = new groovy.text.GStringTemplateEngine()
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

    protected makesafe(TraceRecord record){
        def safefields = [:]
        final stringfields = [
            "status", "hash", "name", "process", "tag", "container", "script", "scratch", "workdir"
        ]
        final longfields = [
            "task_id", "exit", "submit", "start", "attempt", "complete", "duration", "realtime",
            "%cpu", "%mem", "vmem", "rss", "peak_vmem", "peak_rss", "rchar", "wchar",
            "syscr", "syscw", "read_bytes", "write_bytes", "native_id"
        ]
        stringfields.collect{ k ->
            safefields[k] = StringEscapeUtils.escapeJavaScript( record.get(k) as String )
        }
        longfields.collect{ k ->
            safefields[k] = record.get(k) as Long
        }
        safefields
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
}
