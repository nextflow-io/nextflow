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
 * Render pipeline report processes execution
 *
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

    /*
     * Holds an unique index color for each process group
     */
    // final private Map<String,Integer> colorIndexes = new ConcurrentHashMap<>()

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
    }

    /**
     * Save the pending processes and close the trace file
     */
    @Override
    void onFlowComplete() {
        log.debug "Flow completing -- rendering html report"
        endMillis = System.currentTimeMillis()
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


    // final private String REPLACE_STR = '/*REPLACE_WITH_REPORT_DATA*/'

    protected void renderHtml() {

        // make records safe for JS
        def records_safe = [:]
        records.values().eachWithIndex { TraceRecord it, index ->
            records_safe[index] = makesafe(it)
        }

        // render HTML report template
        final tpl_fields = [
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

    // protected StringBuilder renderData() {
    //     final result = new StringBuilder()
    //     result << 'window.data = {\n'
    //     result << '  "trace": [\n'
    //     records.values().eachWithIndex { TraceRecord it, index ->
    //         if( index ) result << ',\n'
    //         append(result, it)
    //     }
    //     result << '\n  ]\n'
    //     result << '};\n'
    //     result
    // }

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

    // protected void append(StringBuilder template, TraceRecord record) {
    //     final name = record.get('name') as String ?: '(unknown)'
    //     final submit = record.get('submit') as Long
    //     final start = record.get('start') as Long
    //     final realtime = record.get('realtime') as Long
    //     final process = record.get('process') as String
    //     final complete = record.get('complete') as Long
    // 
    //     final status = record.get('status') as String
    //     final hash = record.get('hash') as String
    //     final task_id = record.get('task_id') as String
    //     final pct_cpu = record.get('%cpu') as String
    //     final vmem = record.get('vmem') as String
    //     final native_id = record.get('native_id') as String
    //     final exit = record.get('exit') as String
    //     final duration = record.get('duration') as String
    //     final wchar = record.get('wchar') as String
    //     final rchar = record.get('rchar') as String
    //     final rss = record.get('rss') as String
    // 
    //     template << '\n    {\n'
    //     template << "      \"name\": \"${StringEscapeUtils.escapeJavaScript(name)}\", \n"
    //     template << "      \"submit\": \"$submit\", \n"
    //     template << "      \"start\": \"$start\", \n"
    //     template << "      \"realtime\": \"$realtime\", \n"
    //     template << "      \"process\": \"${StringEscapeUtils.escapeJavaScript(process)}\", \n"
    //     template << "      \"complete\": \"$complete\", \n"
    //     template << "      \"status\": \"${StringEscapeUtils.escapeJavaScript(status)}\", \n"
    //     template << "      \"hash\": \"${StringEscapeUtils.escapeJavaScript(hash)}\", \n"
    //     template << "      \"task_id\": \"${StringEscapeUtils.escapeJavaScript(task_id)}\", \n"
    //     template << "      \"pct_cpu\": \"${StringEscapeUtils.escapeJavaScript(pct_cpu)}\", \n"
    //     template << "      \"vmem\": \"${StringEscapeUtils.escapeJavaScript(vmem)}\", \n"
    //     template << "      \"native_id\": \"${StringEscapeUtils.escapeJavaScript(native_id)}\", \n"
    //     template << "      \"exit\": \"${StringEscapeUtils.escapeJavaScript(exit)}\", \n"
    //     template << "      \"duration\": \"${StringEscapeUtils.escapeJavaScript(duration)}\", \n"
    //     template << "      \"wchar\": \"${StringEscapeUtils.escapeJavaScript(wchar)}\", \n"
    //     template << "      \"rchar\": \"${StringEscapeUtils.escapeJavaScript(rchar)}\", \n"
    //     template << "      \"rss\": \"${StringEscapeUtils.escapeJavaScript(rss)}\" \n"
    //     template << '    }'
    // }

    // protected String labelString( TraceRecord record ) {
    //     def result = []
    //     def duration = record.getFmtStr('duration')
    //     def memory = record.getFmtStr('vmem')
    // 
    //     if( duration )
    //         result << duration.toString()
    // 
    //     if( memory )
    //         result <<  memory.toString()
    // 
    //     if( record.cached )
    //         result << 'CACHED'
    // 
    //     result.join(' / ')
    // }

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
