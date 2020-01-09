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

import static nextflow.util.LoggerHelper.*
import static org.fusesource.jansi.Ansi.*

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.util.logging.Slf4j
import jline.TerminalFactory
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
import nextflow.util.Duration
import org.fusesource.jansi.Ansi
import org.fusesource.jansi.AnsiConsole
/**
 * Implements an observer which display workflow
 * execution progress and notifications using
 * ANSI escape codes
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AnsiLogObserver implements TraceObserver {

    static final private String NEWLINE = System.getProperty("line.separator")

    @EqualsAndHashCode
    static class ProcessStats {
        int submitted
        int cached
        int stored
        int failed
        int completed
        boolean terminated
        String hash
        boolean error
        String name
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

    private List<Event> sticky = new ArrayList<>()

    private List<Event> errors = new ArrayList<>()

    private List<Event> warnings = new ArrayList<>()

    private List<Event> infos = new ArrayList<>()

    private Thread renderer

    private Map<String,ProcessStats> processes = new LinkedHashMap()

    private Map<String,Integer> executors = new LinkedHashMap<>()

    private volatile boolean stopped

    private volatile boolean started

    private volatile boolean dirty

    private volatile boolean rendered

    private int printedLines

    private int labelWidth

    private int cols = 1000

    private long startTimestamp

    private long endTimestamp

    private Boolean enableSummary = System.getenv('NXF_ANSI_SUMMARY') as Boolean

    synchronized void appendInfo(String message) {
        boolean warn
        if( isHashLogPrefix(message) && !(warn=message.indexOf('NOTE:')>0) )
            return
        
        if( !started || !processes ) {
            println message
        }
        else if( warn ) {
            warnings << new Event(message)
            dirty = true
        }
        else {
            infos << new Event(message)
            dirty = true
        }
    }

    synchronized void appendWarning(String message) {
        if( !started || !processes )
            printAnsi(message, Color.YELLOW)
        else {
            warnings << new Event(message)
            dirty = true
        }
    }

    synchronized void appendError(String message) {
        if( !started || !processes )
            printAnsi(message, Color.RED)
        else {
            errors << new Event(message)
            dirty = true
            notify()
        }
    }

    synchronized void appendSticky(String message) {
        if( !started || !processes )
            printAnsi(message, Color.RED)
        else {
            sticky << new Event(message)
            dirty = true
            notify()
        }
    }

    protected void render0(dummy) {
        while(!stopped) {
            if( dirty ) renderProgress()
            synchronized (this) {
                wait(200)
            }
        }
        renderProgress()
        renderSummary()
    }

    protected void renderMessages( Ansi term, List<Event> allMessages, Color color=null )  {
        int BLANKS=0
        def itr = allMessages.iterator()
        while( itr.hasNext() ) {
            final event = itr.next()

            // evict old warnings
            final delta = System.currentTimeMillis()-event.timestamp
            if( delta>35_000 ) {
                BLANKS += event.message.count(NEWLINE)+1
                itr.remove()
                continue
            }

            if( color ) {
                term.fg(color).a(event.message).fg(Color.DEFAULT)
            }
            else {
                term.a(event.message)
            }
            term.newline()
        }

        for( int i=0; i<BLANKS; i++ )
            term.newline()
    }

    protected void renderSticky(Ansi term) {
        for( Event ev : sticky ) {
            term.fg(Color.GREEN)
            term.a(Attribute.INTENSITY_BOLD)
            term.a(ev.message)
            term.a(Attribute.INTENSITY_BOLD_OFF)
            term.fg(Color.DEFAULT).newline()
        }
    }

    protected void renderExecutors(Ansi term) {
        int count=0
        def line = ''
        for( Map.Entry<String,Integer> entry : executors ) {
            if( count++ ) line += ","
            line += " $entry.key ($entry.value)"
        }

        if( count ) {
            term.a("executor > " + line)
            term.newline()
        }
    }

    protected void renderProcesses(Ansi term) {
        if( !processes || (!session.isSuccess() && errors && !rendered) ) {
            // prevent to show a useless process progress if there's an error
            // on startup and the execution is terminated
            return
        }

        cols = TerminalFactory.get().getWidth()

        // calc max width
        labelWidth = 0
        for( Map.Entry<String,ProcessStats> entry : processes ) {
            labelWidth = Math.max(labelWidth, entry.value.name.size())
        }

        // render line
        for( Map.Entry<String,ProcessStats> entry : processes ) {
            term.a(line(entry.value))
            term.newline()
        }
        rendered = true
    }

    protected void renderErrors( Ansi term ) {
        for( Event event : errors ) {
            term.fg(Color.RED)
            term.a(event.message)
            term.fg(Color.DEFAULT).newline()
        }
    }

    synchronized protected void renderProgress() {
        dirty = false
        if( printedLines )
            AnsiConsole.out.println ansi().cursorUp(printedLines+1)

        // -- print processes
        final term = ansi()
        renderSticky(term)
        renderExecutors(term)
        renderProcesses(term)
        renderMessages(term, infos)
        renderMessages(term, warnings, Color.YELLOW)
        renderErrors(term)

        final str = term.toString()
        if( str ) {
            printAnsiLines(str)
            printedLines = countNewLines(str)
        }
        else
            printedLines = 0

        AnsiConsole.out.flush()
    }

    protected void renderSummary() {
        final delta = endTimestamp-startTimestamp
        if( enableSummary == false )
            return
        if( enableSummary == null && delta <= 60*1_000 )
            return
        
        WorkflowStats stats = session.isSuccess() ? session.getWorkflowStats() : null
        if( stats && processes.size() ) {
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
            AnsiConsole.out.flush()
        }
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

    protected String fmtWidth(String name, int width, int cols) {
        assert name.size() <= width
        // chop the name string if larger than max cols
        if( name.size() > cols ) {
            return fmtChop(name,cols)
        }
        // otherwise pad with blanks to the expected width
        return name.padRight(Math.min(width,cols))
    }

    protected String fmtChop(String str, int cols) {
        if( str.size() <= cols )
            return str
        return cols>3 ? str[0..(cols-3-1)] + '...' : str[0..cols-1]
    }

    protected String line(ProcessStats stats) {
        final float tot = stats.submitted +stats.cached +stats.stored
        final float com = stats.completed +stats.cached +stats.stored
        final label = fmtWidth(stats.name, labelWidth, Math.max(cols-50, 5))
        final hh = (stats.hash && tot>0 ? stats.hash : '-').padRight(9)

        if( tot == 0  )
            return "[$hh] process > $label -"

        final x = tot ? Math.round(com / tot * 100f) : 0
        final pct = "[${String.valueOf(x).padLeft(3)}%]".toString()

        final numbs = "${(int)com} of ${(int)tot}".toString()
        def result = "[${hh}] process > $label $pct $numbs"
        if( stats.cached )
            result += ", cached: $stats.cached"
        if( stats.stored )
            result += ", stored: $stats.stored"
        if( stats.failed )
            result += ", failed: $stats.failed"
        if( stats.terminated && tot )
            result += stats.error ? ' \u2718' : ' \u2714'
        return fmtChop(result, cols)
    }


    private ProcessStats proc0(String name) {
        def result = processes.get(name)
        if( !result ) {
            result = new ProcessStats()
            result.name = name
            processes.put(name, result)
        }
        return result
    }

    private ProcessStats proc0(TaskHandler handler) {
        proc0(handler.task.processor.name)
    }

    @Override
    void onFlowCreate(Session session){
        this.started = true
        this.session = session
        this.startTimestamp = System.currentTimeMillis()
        AnsiConsole.systemInstall()
        this.renderer = Thread.start('AnsiLogObserver', this.&render0)
    }

    @Override
    void onFlowComplete(){
        stopped = true
        endTimestamp = System.currentTimeMillis()
        renderer.join()
    }

    @Override
    synchronized void onProcessCreate(TaskProcessor processor) {
        proc0(processor.name)
        dirty = true
    }

    /**
     * This method is invoked before a process run is going to be submitted
     * @param handler
     */
    @Override
    synchronized void onProcessSubmit(TaskHandler handler, TraceRecord trace){
        final process = proc0(handler)
        process.submitted++
        process.hash = handler.task.hashLog
        process.name = handler.task.name
        // executor counter
        final exec = handler.task.processor.executor.name
        Integer count = executors[exec] ?: 0
        executors[exec] = count+1
        dirty = true
    }

    @Override
    synchronized void onProcessComplete(TaskHandler handler, TraceRecord trace){
        final process = proc0(handler)
        process.completed++
        process.hash = handler.task.hashLog
        process.name = handler.task.name
        if( handler.task.aborted || handler.task.failed ) {
            process.failed++
            process.error |= !handler.task.errorAction?.soft
        }
        dirty = true
    }

    @Override
    synchronized void onProcessTerminate(TaskProcessor processor) {
        final process = processes.get(processor.name)
        if( process ) {
            process.terminated = true
            dirty = true
        }
    }

    @Override
    synchronized void onProcessCached(TaskHandler handler, TraceRecord trace){
        final process = proc0(handler)
        process.name = handler.task.name
        if( trace ) {
            process.cached++
            process.hash = handler.task.hashLog
        }
        else {
            process.stored++
            process.hash = 'skipped'
        }
        dirty = true
    }

    protected int countNewLines(String str) {
        int result=0
        int gap=0
        while( true ) {
            int p=str.indexOf(NEWLINE)
            if( p==-1 ) {
                result += (str.length()-1).intdiv(cols) +gap
                break
            }
            else {
                result += (p-1).intdiv(cols) +gap
                str = str.substring(p+1)
                if( gap==0 ) gap = 1
            }
        }
        return result
    }
}
