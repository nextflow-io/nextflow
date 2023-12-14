/*
 * Copyright 2013-2023, Seqera Labs
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
import jline.TerminalFactory
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.util.Duration
import nextflow.util.Threads
import org.fusesource.jansi.Ansi
import org.fusesource.jansi.AnsiConsole
import static nextflow.util.LoggerHelper.isHashLogPrefix
import static org.fusesource.jansi.Ansi.Attribute
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

    static final private String NEWLINE = '\n'

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

    private Map<String,Integer> executors = new LinkedHashMap<>()

    private volatile boolean stopped

    private volatile boolean started

    private volatile long changeTimestamp

    private volatile long lastRendered

    private volatile boolean rendered

    private int printedLines

    private int gapLines

    private int labelWidth

    private volatile int cols = 80

    private volatile int rows = 24

    private long startTimestamp

    private long endTimestamp

    private long lastWidthReset

    private Boolean enableSummary = System.getenv('NXF_ANSI_SUMMARY') as Boolean

    private final int WARN_MESSAGE_TIMEOUT = 35_000

    private WorkflowStatsObserver statsObserver

    private void markModified() {
        changeTimestamp = System.currentTimeMillis()
    }

    boolean getStarted() { started }

    boolean  getStopped() { stopped }

    private boolean hasProgressChanges() {
        final long progress = statsObserver.changeTimestamp ?: 0
        final last = progress ? Math.max(progress,changeTimestamp) : changeTimestamp
        if( last != lastRendered ) {
            lastRendered = last
            return true
        }
        return false
    }

    synchronized void appendInfo(String message) {
        if( message==null )
            return
        boolean warn
        if( isHashLogPrefix(message) && !(warn=message.indexOf('NOTE:')>0) )
            return
        
        if( !started || !statsObserver.hasProgressRecords() ) {
            println message
        }
        else if( warn ) {
            warnings << new Event(message)
            markModified()
        }
        else {
            infos << new Event(message)
            markModified()
        }
    }

    synchronized void appendWarning(String message) {
        if( message==null )
            return
        if( !started || !statsObserver.hasProgressRecords() )
            printAnsi(message, Color.YELLOW)
        else {
            warnings << new Event(message)
            markModified()
        }
    }

    synchronized void appendError(String message) {
        if( message==null )
            return
        if( !started || !statsObserver.hasProgressRecords() ) {
            printAnsi(message, Color.RED)
        }
        else {
            errors << new Event(message)
            markModified()
            notify()
        }
    }

    synchronized void appendSticky(String message) {
        if( message==null )
            return
        if( !started || !statsObserver.hasProgressRecords() )
            printAnsi(message, Color.GREEN)
        else {
            sticky << new Event(message)
            markModified()
            notify()
        }
    }

    protected void render0(dummy) {
        while(!stopped) {
            if( hasProgressChanges() )
                renderProgress(statsObserver.quickStats)
            synchronized (this) {
                wait(200)
            }
        }
        // 
        final stats = statsObserver.getStats()
        renderProgress(stats)
        renderSummary(stats)
    }

    protected void renderMessages( Ansi term, List<Event> allMessages, Color color=null )  {
        int BLANKS=0
        def itr = allMessages.iterator()
        while( itr.hasNext() ) {
            final event = itr.next()

            // evict old warnings
            final delta = System.currentTimeMillis()-event.timestamp
            if( delta>WARN_MESSAGE_TIMEOUT ) {
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

    protected String getExecutorName(String key) {
        session.getExecutorFactory().getDisplayName(key)
    }
    
    protected void renderExecutors(Ansi term) {
        int count=0
        def line = ''
        for( Map.Entry<String,Integer> entry : executors ) {
            if( count++ ) line += ","
            line += " ${getExecutorName(entry.key)} ($entry.value)"
        }

        if( count ) {
            term.a("executor > " + line)
            term.newline()
        }
    }

    protected void renderProcesses(Ansi term, WorkflowStats stats) {
        def processes = stats.getProcesses()
        if( !processes || (!session.isSuccess() && errors && !rendered) ) {
            // prevent to show a useless process progress if there's an error
            // on startup and the execution is terminated
            return
        }

        cols = TerminalFactory.get().getWidth()
        rows = TerminalFactory.get().getHeight()

        // calc max width
        final now = System.currentTimeMillis()
        if( now-lastWidthReset>20_000 )
            labelWidth = 0

        final lastWidth = labelWidth
        for( ProgressRecord entry : processes ) {
            labelWidth = Math.max(labelWidth, entry.taskName.size())
        }
        if( lastWidth != labelWidth )
            lastWidthReset = now

        // render line
        def renderedLines = 0
        def skippedLines = 0
        for( ProgressRecord entry : processes ) {
            if( renderedLines <= rows - 5 || entry.getTotalCount() > 0 ) {
                term = line(entry, term)
                term.newline()
                renderedLines += 1
            }
            else {
                skippedLines += 1
            }
        }
        if( skippedLines > 0 )
            term.a("Plus $skippedLines more processes waiting for tasks…").newline()
        rendered = true
    }

    protected void renderErrors( Ansi term ) {
        for( Event event : errors ) {
            term.fg(Color.RED)
            term.a(event.message)
            term.fg(Color.DEFAULT).newline()
        }
    }

    synchronized protected void renderProgress(WorkflowStats stats) {
        if( printedLines )
            AnsiConsole.out().println ansi().cursorUp(printedLines+gapLines+1)

        // -- print processes
        final term = ansi()
        renderSticky(term)
        renderExecutors(term)
        renderProcesses(term, stats)
        renderMessages(term, infos)
        renderMessages(term, warnings, Color.YELLOW)
        renderErrors(term)

        final str = term.toString()
        final count = printAndCountLines(str)
        AnsiConsole.out().flush()

        // usually the gap should be negative because `count` should be greater or equal
        // than the previous `printedLines` value (the output should become longer)
        // otherwise cleanup the remaining lines
        gapLines = printedLines > count ? printedLines-count : 0
        if( gapLines>0 ) for(int i=0; i<gapLines; i++ )
            AnsiConsole.out().print(ansi().eraseLine().newline())
        // at the end update the value of printed lines
        printedLines = count
    }

    protected int printAndCountLines(String str) {
        if( str ) {
            printAnsiLines(str)
            return str.count(NEWLINE)
        }
        else
            return 0
    }

    protected void renderSummary(WorkflowStats stats) {
        final delta = endTimestamp-startTimestamp
        if( enableSummary == false )
            return
        if( enableSummary == null && delta <= 60*1_000 )
            return
        
        if( session.isSuccess() && stats.progressLength>0 ) {
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
            AnsiConsole.out().flush()
        }
    }

    protected void printAnsi(String message, Color color=null, boolean bold=false) {
        def fmt = ansi()
        if( color ) fmt = fmt.fg(color)
        if( bold ) fmt = fmt.bold()
        fmt = fmt.a(message)
        if( bold ) fmt = fmt.boldOff()
        if( color ) fmt = fmt.fg(Color.DEFAULT)
        AnsiConsole.out().println(fmt.eraseLine())
    }
    
    protected void printAnsiLines(String lines) {
        final text = lines
                .replace('\r','')
                .replace(NEWLINE, ansi().eraseLine().toString() + NEWLINE)
        AnsiConsole.out().print(text)
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
        return str.take(3) + '…' + str.takeRight(cols-1-3)
    }

    protected Ansi line(ProgressRecord stats, Ansi term) {
        final float tot = stats.getTotalCount()
        final float com = stats.getCompletedCount()
        final label = fmtWidth(stats.taskName, labelWidth, Math.max(cols-50, 5))
        final tagMatch = label =~ /( \(.+\) *)$/
        final labelTag = tagMatch ? tagMatch.group(1) : ''
        final labelNoTag = label.replaceFirst(/ \(.+\) *$/, "")
        final hh = (stats.hash && tot>0 ? stats.hash : '-').padRight(9)

        final x = tot ? Math.floor(com / tot * 100f).toInteger() : 0
        final pct = "${String.valueOf(x).padLeft(3)}%".toString()
        final numbs = " ${(int)com} of ${(int)tot}".toString()

        // Hash: []
        term = term.fg(Color.BLACK).a('[').fg(Color.BLUE).a(hh).fg(Color.BLACK).a('] ')

        // Label: process > sayHello
        if( cols > 180 )
            term = term.a('process > ')
        term = term.reset().a(labelNoTag).fg(Color.YELLOW).a(labelTag).reset()

        // No tasks
        if( tot == 0 ) {
            term = term.a(' -')
            return term
        }

        // Progress: [  0%] 0 of 10
        if( cols > 120 ) {
            term = term
                .fg(Color.BLACK)
                .a(' [')
                .fg(pct == '100%' ? Color.GREEN : Color.YELLOW)
                .a(pct)
                .fg(Color.BLACK)
                .a(']')
                .reset()
        }
        else {
            term = term.fg(Color.BLACK).a(' |').reset()
        }
        term = term.a(numbs)

        // Number of tasks:
        if( stats.cached )
            term = term.fg(Color.BLACK).a(", cached: $stats.cached").reset()
        if( stats.stored )
            term = term.a(", stored: $stats.stored")
        if( stats.failed )
            term = term.a(", failed: $stats.failed")
        if( stats.retries )
            term = term.a(", retries: $stats.retries")
        if( stats.terminated && tot ) {
            term = stats.errored
                ? term.fg(Color.RED).a(' \u2718' ).reset()
                : term.fg(Color.GREEN).a(' \u2714 ').reset()
        }
        return term
    }

    @Override
    void onFlowCreate(Session session){
        this.started = true
        this.session = session
        this.statsObserver = session.statsObserver
        this.startTimestamp = System.currentTimeMillis()
        AnsiConsole.systemInstall()
        this.renderer = Threads.start('AnsiLogObserver', this.&render0)
    }

    @Override
    void onFlowComplete(){
        stopped = true
        endTimestamp = System.currentTimeMillis()
        renderer.join()
    }

    /**
     * This method is invoked before a process run is going to be submitted
     * @param handler
     */
    @Override
    synchronized void onProcessSubmit(TaskHandler handler, TraceRecord trace){
        // executor counter
        final exec = handler.task.processor.executor.name
        Integer count = executors[exec] ?: 0
        executors[exec] = count+1
        markModified()
    }

    void forceTermination() {
        stopped = true
        endTimestamp = System.currentTimeMillis()
    }
}
