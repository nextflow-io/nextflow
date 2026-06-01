/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cli

import java.nio.charset.StandardCharsets
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.util.regex.Pattern

import ch.qos.logback.classic.Level
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.SysEnv
import nextflow.cli.LogFileFormatter.EntryHeader
import nextflow.exception.AbortOperationException
import nextflow.util.HistoryFile

/**
 * Implements the `logfile` command to print the contents of a Nextflow log file,
 * resolved either by run name (via {@link HistoryFile}) or by an explicit file path.
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Print the contents of a Nextflow log file")
class CmdLogFile extends CmdBase {

    static final public NAME = 'logfile'

    static final String LOG_FILE_NAME = '.nextflow.log'

    static final String DEFAULT_QUERY = 'last'

    static final List<String> LEVEL_NAMES = ['TRACE', 'DEBUG', 'INFO', 'WARN', 'ERROR']

    /** Filename suffixes scanned when resolving a run name to a (possibly rotated) log file. */
    static final List<String> LOG_FILE_SUFFIXES = ['', '.1', '.2', '.3', '.4', '.5', '.6', '.7', '.8', '.9']

    /**
     * Cap on lines scanned when locating a session UUID. The UUID is logged near the top of the
     * file, but verbose plugin/classpath logging at -trace/-debug can push it down a few hundred
     * lines, so the cap is generous to avoid spurious "cannot find log" resolutions.
     */
    static final int SESSION_ID_SCAN_LIMIT = 1000

    static final long FOLLOW_POLL_MS = 250

    /** Matches an ANSI CSI escape sequence (color and style codes) found in input. */
    static final Pattern ANSI_CSI = ~/\[[\d;?]*[@-~]/

    @Parameter(description = 'Run name, session id, or path to a .nextflow.log file')
    List<String> args

    @Parameter(names = ['-l','-level'], description = 'Minimum log level to print: TRACE, DEBUG, INFO, WARN, ERROR (no filtering by default)')
    String level

    @Parameter(names = ['-f','-follow'], description = 'Output appended data as the file grows', arity = 0)
    boolean follow

    @Parameter(names = ['-n','-lines'], description = 'Print only the last N log entries')
    Integer lines

    @Parameter(names = ['-no-ansi'], description = 'Do not add color/style codes to the output (also disabled when NO_COLOR is set or stdout is not a terminal)', arity = 0)
    boolean noAnsi

    @Parameter(names = ['-keep-ansi'], description = 'Preserve ANSI escape codes found inside log content (default: strip them)', arity = 0)
    boolean keepAnsi

    @Parameter(names = ['-no-pager'], description = 'Do not pipe the output through a pager (default: pipe through $NXF_PAGER, then $PAGER, then `less -FR` when stdout is a terminal)', arity = 0)
    boolean noPager

    @PackageScope HistoryFile history

    @PackageScope Path currentDir

    /** Test seam: when non-null, overrides TTY detection during color resolution. */
    @PackageScope Boolean colorOverride

    private PrintStream out = System.out

    private boolean stripAnsi
    private LogFileFormatter formatter

    @Override
    final String getName() { NAME }

    @Override
    void run() {
        if( args && args.size() > 1 )
            throw new AbortOperationException("Too many arguments — expected a single run name or file path")
        if( lines != null && lines < 0 )
            throw new AbortOperationException("Option `-n` must be a non-negative integer")

        stripAnsi = !keepAnsi
        final usePager = shouldUsePager()
        formatter = new LogFileFormatter(resolveUseColor(usePager))

        final threshold = parseLevel(level)
        final path = resolveLogFile(args ? args[0] : null)

        if( follow ) {
            // Snapshot the size once so the tail read and the follow seek share the same byte
            // boundary — otherwise lines appended between them would be printed by neither.
            final long startPos = Files.size(path)
            if( lines != null )
                printLastN(path, threshold, lines, startPos)
            followFile(path, threshold, startPos)
            return
        }

        if( usePager )
            runWithPager { doPrint(path, threshold) }
        else
            doPrint(path, threshold)
    }

    private void doPrint(Path path, Level threshold) {
        if( lines != null )
            printLastN(path, threshold, lines)
        else
            printAll(path, threshold)
    }

    private boolean resolveUseColor(boolean usePager) {
        if( noAnsi )
            return false
        // Honor the same global gates as ColorUtil.isAnsiEnabled(): NO_COLOR disables, then
        // NXF_ANSI_LOG forces on/off. These take precedence over pager/TTY detection so a user
        // who disables color globally also gets uncolored `logfile` output.
        if( SysEnv.get('NO_COLOR') != null )
            return false
        final ansiLog = SysEnv.get('NXF_ANSI_LOG')?.trim()
        if( ansiLog )
            return Boolean.parseBoolean(ansiLog)
        if( colorOverride != null )
            return colorOverride
        if( usePager )
            return true  // pager handles the real terminal — keep colors
        return System.console() != null
    }

    @PackageScope
    boolean shouldUsePager() {
        if( noPager ) return false
        if( follow )  return false   // pagers buffer; tail-follow needs an unbuffered stream
        return System.console() != null
    }

    @PackageScope
    static List<String> resolvePagerCommand() {
        final nxfPager = SysEnv.get('NXF_PAGER')?.trim()
        if( nxfPager ) return parsePagerCommand(nxfPager)
        final pager = SysEnv.get('PAGER')?.trim()
        if( pager ) return parsePagerCommand(pager)
        return ['less', '-FR']
    }

    @PackageScope
    static List<String> parsePagerCommand(String value) {
        final parts = value.trim().tokenize() as List<String>
        if( parts.isEmpty() )
            return null
        // Bare `less` (no args, with or without a path prefix) gets `-FR` appended so colors
        // render (-R) and short output exits without entering the pager (-F). We omit `-X`
        // (which git's default includes) because it disables the alternate screen buffer that
        // most terminals rely on for mouse-wheel-to-arrow-key translation.
        final cmd = parts[0]
        final base = cmd.contains('/') ? cmd.substring(cmd.lastIndexOf('/') + 1) : cmd
        if( parts.size() == 1 && base == 'less' )
            parts.add('-FR')
        return parts
    }

    private void runWithPager(Runnable printAction) {
        final cmd = resolvePagerCommand()
        if( !cmd ) {
            printAction.run()
            return
        }
        Process proc
        try {
            proc = new ProcessBuilder(cmd)
                .redirectOutput(ProcessBuilder.Redirect.INHERIT)
                .redirectError(ProcessBuilder.Redirect.INHERIT)
                .start()
        }
        catch( IOException e ) {
            log.debug "Could not start pager '${cmd.join(' ')}': ${e.message}"
            printAction.run()
            return
        }
        final originalOut = this.out
        // Buffer the pipe to the pager so each autoflush'd println drains a coalesced
        // buffer rather than emitting multiple small writes per line.
        this.out = new PrintStream(new BufferedOutputStream(proc.outputStream), true, StandardCharsets.UTF_8)
        try {
            printAction.run()
        }
        finally {
            try { this.out.close() } catch( Exception ignored ) {}
            try { proc.waitFor() } catch( InterruptedException ignored) {}
            this.out = originalOut
        }
    }

    private Path resolveLogFile(String arg) {
        if( arg ) {
            final candidate = Paths.get(arg)
            if( Files.exists(candidate) ) {
                if( !Files.isRegularFile(candidate) )
                    throw new AbortOperationException("Not a regular file: $arg")
                return candidate
            }
        }

        final query = arg?.trim() ?: DEFAULT_QUERY
        final hist = history ?: HistoryFile.DEFAULT
        if( !hist.exists() || hist.empty() )
            throw new AbortOperationException("It looks like no pipeline was executed (execution history is empty)")

        final records = hist.findByIdOrName(query)
        if( !records )
            throw new AbortOperationException("Unknown run name or session id: $query")

        final record = records.last()
        if( record.sessionId == null )
            throw new AbortOperationException("History record for run '${record.runName ?: query}' has no session id — pass the log file path explicitly")
        final sessionId = record.sessionId.toString()
        final dir = currentDir ?: Paths.get('.')
        for( String suffix : LOG_FILE_SUFFIXES ) {
            final candidate = dir.resolve("${LOG_FILE_NAME}${suffix}".toString())
            if( Files.exists(candidate) && fileContainsSessionId(candidate, sessionId) )
                return candidate
        }

        throw new AbortOperationException(
            "Cannot find a ${LOG_FILE_NAME} file in the current directory for run '${record.runName ?: query}' (session ${sessionId}). " +
            "Pass the log file path explicitly.")
    }

    private static boolean fileContainsSessionId(Path file, String sessionId) {
        try {
            return Files.newBufferedReader(file, StandardCharsets.UTF_8).withCloseable { reader ->
                String line
                int count = 0
                while( count++ < SESSION_ID_SCAN_LIMIT && (line = reader.readLine()) != null ) {
                    if( line.contains(sessionId) )
                        return true
                }
                return false
            }
        }
        catch( IOException e ) {
            log.debug "Could not scan file ${file}: ${e.message}"
            return false
        }
    }

    @PackageScope
    static Level parseLevel(String value) {
        if( !value )
            return null
        final str = value.trim().toUpperCase()
        // Level.toLevel accepts extras like ALL/OFF — gate on our supported set first
        if( !(str in LEVEL_NAMES) )
            throw new AbortOperationException("Invalid log level: '${value}' — expected one of ${LEVEL_NAMES.join(', ')}")
        return Level.toLevel(str)
    }

    private static boolean passesFilter(Level entryLevel, Level threshold) {
        if( threshold == null )
            return true
        // treat unknown (continuation before any entry-start) as INFO so it isn't silently dropped at default thresholds
        return (entryLevel ?: Level.INFO).isGreaterOrEqual(threshold)
    }

    private String stripContentAnsi(String line) {
        stripAnsi ? ANSI_CSI.matcher(line).replaceAll('') : line
    }

    private void printAll(Path path, Level threshold) {
        boolean keep = passesFilter(null, threshold)
        Level currentLevel = null
        Files.newBufferedReader(path, StandardCharsets.UTF_8).withCloseable { reader ->
            String line
            while( (line = reader.readLine()) != null ) {
                line = stripContentAnsi(line)
                final header = LogFileFormatter.parseEntryHeader(line)
                if( header != null ) {
                    currentLevel = header.level
                    keep = passesFilter(currentLevel, threshold)
                }
                if( keep ) {
                    out.println(formatter.format(line, currentLevel, header))
                    if( out.checkError() ) break   // pager exited (broken pipe) — stop reading
                }
            }
        }
    }

    private void printLastN(Path path, Level threshold, int n, long endPos = Long.MAX_VALUE) {
        if( n == 0 )
            return
        final ArrayDeque<String> buffer = new ArrayDeque<>(n)
        StringBuilder pending = null
        Level pendingLevel = null
        boolean pendingKeep = passesFilter(null, threshold)

        final flush = {
            if( pending != null && pendingKeep ) {
                if( buffer.size() == n )
                    buffer.removeFirst()
                buffer.addLast(pending.toString())
            }
        }

        boundedReader(path, endPos).withCloseable { reader ->
            String line
            while( (line = reader.readLine()) != null ) {
                line = stripContentAnsi(line)
                final header = LogFileFormatter.parseEntryHeader(line)
                if( header != null ) {
                    flush.call()
                    pending = new StringBuilder().append(formatter.format(line, header.level, header))
                    pendingKeep = passesFilter(header.level, threshold)
                    pendingLevel = header.level
                }
                else if( pending == null ) {
                    pending = new StringBuilder().append(formatter.format(line, null, null))
                    pendingLevel = null
                }
                else {
                    pending.append('\n').append(formatter.format(line, pendingLevel, null))
                }
            }
        }
        flush.call()

        for( String entry : buffer ) {
            out.println(entry)
            if( out.checkError() ) break   // pager exited
        }
    }

    /**
     * Open a reader over {@code path}, optionally bounded to the first {@code endPos} bytes so a
     * tail read in follow mode stops exactly where {@link #followFile} resumes.
     */
    private static BufferedReader boundedReader(Path path, long endPos) {
        final InputStream raw = Files.newInputStream(path)
        final InputStream src = endPos == Long.MAX_VALUE ? raw : new BoundedInputStream(raw, endPos)
        return new BufferedReader(new InputStreamReader(src, StandardCharsets.UTF_8))
    }

    /** Caps the number of bytes read from a delegate stream. */
    @CompileStatic
    private static class BoundedInputStream extends InputStream {
        private final InputStream delegate
        private long remaining

        BoundedInputStream(InputStream delegate, long limit) {
            this.delegate = delegate
            this.remaining = limit
        }

        @Override
        int read() throws IOException {
            if( remaining <= 0 ) return -1
            final b = delegate.read()
            if( b >= 0 ) remaining--
            return b
        }

        @Override
        int read(byte[] buf, int off, int len) throws IOException {
            if( remaining <= 0 ) return -1
            final n = delegate.read(buf, off, (int) Math.min((long) len, remaining))
            if( n > 0 ) remaining -= n
            return n
        }

        @Override
        void close() throws IOException { delegate.close() }
    }

    // RandomAccessFile.readLine() reads bytes as ISO-8859-1. Nextflow log content is essentially
    // ASCII (logger names, ISO timestamps, English messages), so this is equivalent to UTF-8 in
    // practice. Non-ASCII characters in messages may be displayed incorrectly in follow mode.
    private void followFile(Path path, Level threshold, long startPos) {
        boolean keep = passesFilter(null, threshold)
        Level currentLevel = null
        RandomAccessFile raf = null
        try {
            raf = new RandomAccessFile(path.toFile(), 'r')
            raf.seek(startPos)
            while( !Thread.currentThread().isInterrupted() ) {
                if( raf.length() < raf.getFilePointer() ) {
                    // file was truncated or rotated — reopen from the start
                    raf.close()
                    raf = new RandomAccessFile(path.toFile(), 'r')
                    keep = passesFilter(null, threshold)
                    currentLevel = null
                }
                String line
                boolean read = false
                while( (line = raf.readLine()) != null ) {
                    read = true
                    line = stripContentAnsi(line)
                    final header = LogFileFormatter.parseEntryHeader(line)
                    if( header != null ) {
                        currentLevel = header.level
                        keep = passesFilter(currentLevel, threshold)
                    }
                    if( keep )
                        out.println(formatter.format(line, currentLevel, header))
                }
                if( read )
                    out.flush()
                Thread.sleep(FOLLOW_POLL_MS)
            }
        }
        catch( InterruptedException e ) {
            Thread.currentThread().interrupt()
        }
        finally {
            try { raf?.close() } catch( Exception ignored ) {}
        }
    }
}
