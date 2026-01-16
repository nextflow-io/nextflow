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
import java.util.concurrent.ConcurrentHashMap

import groovy.text.GStringTemplateEngine
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.processor.TaskId
import nextflow.trace.config.TimelineConfig
import nextflow.trace.event.TaskEvent
import nextflow.util.Duration
import nextflow.util.TestOnly
import org.apache.commons.lang3.StringEscapeUtils
/**
 * Render pipeline timeline processes execution
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TimelineObserver implements TraceObserverV2 {

    /**
     * Holds the the start time for tasks started/submitted but not yet completed
     */
    final private Map<TaskId,TraceRecord> records = new LinkedHashMap<>()

    /*
     * Holds an unique index color for each process group
     */
    final private Map<String,Integer> colorIndexes = new ConcurrentHashMap<>()

    private Path reportFile

    private long beginMillis

    private long startMillis

    private long endMillis

    private boolean overwrite

    TimelineObserver(TimelineConfig config) {
        this.reportFile = FileHelper.asPath(config.file)
        this.overwrite = config.overwrite
    }

    @TestOnly
    protected TimelineObserver() {}

    /**
     * Create the trace file, in file already existing with the same name it is
     * "rolled" to a new file
     */
    @Override
    void onFlowCreate(Session session) {
        beginMillis = startMillis = System.currentTimeMillis()
        if( Files.exists(reportFile) && !overwrite )
            throw new AbortOperationException("Timeline file already exists: ${reportFile.toUriString()} -- enable the 'timeline.overwrite' option in your config file to overwrite existing files")
    }

    /**
     * Save the pending processes and close the trace file
     */
    @Override
    void onFlowComplete() {
        log.debug "Workflow completed -- rendering execution timeline"
        endMillis = System.currentTimeMillis()
        try {
            renderHtml()
        }
        catch (Exception e) {
            log.warn "Failed to render execution timeline -- see the log file for details", e
        }
    }

    @Override
    void onTaskSubmit(TaskEvent event) {
        synchronized (records) {
            records[ event.trace.taskId ] = event.trace
        }
    }

    @Override
    void onTaskStart(TaskEvent event) {
        synchronized (records) {
            records[ event.trace.taskId ] = event.trace
        }
    }

    @Override
    void onTaskComplete(TaskEvent event) {
        final taskId = event.handler.task.id
        if( !event.trace ) {
            log.debug "Profile warn: Unable to find record for task_run with id: ${taskId}"
            return
        }

        // remove the record from the current records
        synchronized (records) {
            records[ event.trace.taskId ] = event.trace
        }
    }

    @Override
    void onTaskCached(TaskEvent event) {

        // event was triggered by a stored task, ignore it
        if( !event.trace ) {
            return
        }

        // remove the record from the current records
        synchronized (records) {
            records[ event.trace.taskId ] = event.trace
        }

        beginMillis = Math.min( beginMillis, event.trace.get('submit') as long )
    }

    protected void renderHtml() {
        // render HTML timeline template
        final tpl_fields = [
                payload : renderData(),
                assets_js : [
                        readTemplate('assets/jquery-3.2.1.min.js'),
                        readTemplate('assets/d3.v3.min.js'),
                        readTemplate('assets/d3-timeline.min.js'),
                        readTemplate('assets/TimelineTemplate.js')
                ]
        ]

        final tpl = readTemplate('TimelineTemplate.html')
        final engine = new GStringTemplateEngine()
        final html_template = engine.createTemplate(tpl)
        final html_output = html_template.make(tpl_fields).toString()

        // make sure the parent path exists
        final parent = reportFile.getParent()
        if( parent )
            Files.createDirectories(parent)

        final writer = TraceHelper.newFileWriter(reportFile, overwrite, 'Timeline')
        writer.withWriter { w -> w << html_output }
        writer.close()
    }

    protected String renderData() {
        // Returns a JSON formatted string
        final result = new StringBuilder()
        final indent = "    ";
        result << '{\n'
        result << indent << '"elapsed": "' << new Duration(endMillis-startMillis).toString() << '",\n'
        result << indent << '"beginningMillis": ' << beginMillis << ',\n'
        result << indent << '"endingMillis": ' << endMillis << ',\n'
        result << indent << '"processes": [\n'
        records.values().eachWithIndex { TraceRecord it, index ->
            if( index ) result << ',\n'
            append(result, it)
        }
        result << '\n' << indent << ']\n'
        result << '}\n'
        return result.toString()
    }

    protected void append(StringBuilder template, TraceRecord record) {
        final name = record.get('name') as String ?: '(unknown)'
        final submit = record.get('submit') as Long
        final start = record.get('start') as Long
        final realtime = record.get('realtime') as Long
        final process = record.get('process') as String
        final complete = record.get('complete') as Long
        final index = colorIndexes.getOrCreate(process) { colorIndexes.size() }
        final indent = "    ";

        template << indent << indent << '{'
        template << "\"label\": \"${StringEscapeUtils.escapeEcmaScript(name)}\", "
        template << "\"cached\": ${record.cached}, "
        template << "\"index\": $index, "
        template << "\"times\": ["

        if( submit && start ) {
            template << "{\"starting_time\": $submit, \"ending_time\": $start}"

            if( start && realtime ) {
                final label = StringEscapeUtils.escapeEcmaScript(labelString(record))
                final ending = start+realtime
                template << ", {\"starting_time\": $start, \"ending_time\": $ending, \"label\": \"$label\"}"

                if( complete && ending < complete ) {
                    template << ", {\"starting_time\": $ending, \"ending_time\": $complete}"
                }
            }
        }

        template << "]}"
    }

    protected String labelString( TraceRecord record ) {
        final result = []
        final duration = record.getFmtStr('duration')
        final mem = record.getFmtStr('peak_rss')

        if( duration )
            result << duration.toString()

        if( mem )
            result << mem.toString()

        if( record.cached )
            result << 'CACHED'

        result.join(' / ')
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

    @Override
    boolean enableMetrics() {
        return true
    }
}
