/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.util.Duration
import org.fusesource.jansi.Ansi
import org.fusesource.jansi.AnsiConsole
import static org.fusesource.jansi.Ansi.Color
import static org.fusesource.jansi.Ansi.ansi

/**
 * Implements an observer which display workflow
 * execution progress and notifications using
 * ANSI escape codes
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AnsiLogObserver implements TraceObserver {

    static final private String NEWLINE = System.getProperty("line.separator")

    static class ProcessStats {
        int submitted
        int cached
        int failed
        int completed
        String hash
        boolean error
        boolean changed
    }

    static class Event {
        String message
        long timestamp

        Event(String message) {
            this.message = message
            this.timestamp = System.currentTimeMillis()
        }

        @Override
        String toString() { return message }
    }

    private Session session

    private List<Event> errors = new ArrayList<>()

    private List<Event> warnings = new ArrayList<>()

    private List<Event> infos = new ArrayList<>()

    private Thread renderer

    private Map<String,ProcessStats> processes = new LinkedHashMap()

    private volatile boolean stopped

    private volatile boolean started

    private int printedLines

    private int maxNameLength

    private WorkflowStats stats

    private long startTimestamp

    private long endTimestamp

    synchronized void appendInfo(String message) {
        if( message.startsWith('[') )
            return
        
        if( !started || !processes ) {
            println message
        }
        else {
            if( message.indexOf('NOTE:')!=-1 )
                warnings << new Event(message)
            else
                infos << new Event(message)
        }
    }

    synchronized void appendWarning(String message) {
        if( !started || !processes )
            printAnsi(message, Color.YELLOW)
        else
            warnings << new Event(message)
    }

    synchronized void appendError(String message) {
        if( !started || !processes )
            printAnsi(message, Color.RED)
        else
            errors << new Event(message)
    }

    protected void render0(dummy) {
        while(!stopped) {
            renderProgress()
            sleep 150
        }
    }

    protected void renderMessages( Ansi term, List<Event> allMessages, Color color=null )  {
        def itr = allMessages.iterator()
        while( itr.hasNext() ) {
            final event = itr.next()

            if( color ) {
                term.fg(color).a(event.message).fg(Color.DEFAULT)
            }
            else {
                term.a(event.message)
            }
            term.newline()

            // evict old warnings
            final delta = System.currentTimeMillis()-event.timestamp
            if( delta>60_000 )
                itr.remove()
        }
    }

    protected void renderProcesses(Ansi term) {
        for( Map.Entry<String,ProcessStats> entry : processes ) {
            term.a(line(entry.key,entry.value))
            term.newline()
        }
    }

    protected void renderErrors( Ansi term ) {
        for( Event event : errors ) {
            term.fg(Color.RED)
            term.a(event.message)
            term.fg(Color.DEFAULT).newline()
        }
    }

    synchronized protected void renderProgress() {
        if( printedLines )
            AnsiConsole.out.println ansi().cursorUp(printedLines+1)

        // -- print processes
        final term = ansi()
        renderProcesses(term)
        renderMessages(term, infos)
        renderMessages(term, warnings, Color.YELLOW)
        renderErrors(term)

        final str = term.toString()
        if( str ) {
            printAnsiLines(str)
            printedLines = str.count(NEWLINE)
        }
        else
            printedLines = 0

        if( stats ) {
            final delta = endTimestamp-startTimestamp
            def report = ""
            report += "Completed at: ${new Date(endTimestamp).format('dd-MMM-yyyy HH:mm:ss')}\n"
            report += "Duration    : ${new Duration(delta)}\n"
            report += "CPU hours   : ${stats.getComputeTimeFmt()}\n"
            report += "Succeeded   : ${stats.succeedCountFmt}\n"
            if( stats.cachedCount )
                report += "Cached      : ${stats.cachedCountFmt}\n"
            if( stats.ignoredCount )
                report += "Ignored     : ${stats.ignoredCountFmt}\n"
            if( stats.failedCount )
                report += "Failed      : ${stats.failedCountFmt}\n"

            printAnsi(report, Color.GREEN, true)
        }

        AnsiConsole.out.flush()
    }

    protected void printAnsi(String message, Ansi.Color color=null, boolean bold=false) {
        def fmt = ansi()
        if( color ) fmt = fmt.fg(color)
        if( bold ) fmt = fmt.bold()
        fmt = fmt.a(message)
        if( bold ) fmt = fmt.boldOff()
        if( color ) fmt = fmt.fg(Color.DEFAULT)
        AnsiConsole.out.println(fmt.eraseLine())
    }
    
    protected void printAnsiLines(String lines) {
        final text = lines.replace(NEWLINE,  ansi().eraseLine().toString() + NEWLINE)
        AnsiConsole.out.print(text)
    } 

    protected String line(String name, ProcessStats stats) {
        final int tot = stats.submitted + stats.cached
        final int com = stats.completed + stats.cached
        def x = Math.round((float)(com / tot) * 100)
        final pct = "[${String.valueOf(x).padLeft(3)}%]".toString()
        final label = name.padRight(maxNameLength)
        final numbs = "${com} of ${tot}".toString()
        def result = "[$stats.hash] $label: $pct $numbs"
        if( stats.cached )
            result += ", cached: $stats.cached"
        if( stats.failed )
            result += ", failed: $stats.failed"
        return result
    }


    private ProcessStats p0(String name) {
        def result = processes.get(name)
        if( !result ) {
            result = new ProcessStats()
            processes.put(name, result)
            maxNameLength = Math.max(maxNameLength, name.size())
        }
        return result
    }

    private ProcessStats p0(TaskHandler handler) {
        p0(handler.task.processor.name)
    }

    @Override
    void onFlowStart(Session session){
        this.started = true
        this.session = session
        this.startTimestamp = System.currentTimeMillis()
        AnsiConsole.systemInstall()
        this.renderer = Thread.start('AnsiLogObserver', this.&render0)
    }

    @Override
    void onFlowComplete(){
        this.stopped = true
        this.stats = session.isSuccess() ? session.getWorkflowStats() : null
        this.endTimestamp = System.currentTimeMillis()
        renderProgress()
    }

    /**
     * This method is invoked before a process run is going to be submitted
     * @param handler
     */
    @Override
    synchronized void onProcessSubmit(TaskHandler handler, TraceRecord trace){
        final process = p0(handler)
        process.submitted++
        process.hash = handler.task.hashLog
        process.changed = true
    }

    @Override
    synchronized void onProcessComplete(TaskHandler handler, TraceRecord trace){
        final process = p0(handler)
        process.completed++
        process.hash = handler.task.hashLog
        process.changed = true
        if( handler.getStatusString() != 'COMPLETED' ) {
            process.failed++
            process.error = true
        }
    }

    @Override
    synchronized void onProcessCached(TaskHandler handler, TraceRecord trace){
        final process = p0(handler)
        process.cached++
        process.hash = handler.task.hashLog
        process.changed = true
    }

}
