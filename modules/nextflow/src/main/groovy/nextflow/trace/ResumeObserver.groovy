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

package nextflow.trace

import java.nio.file.Paths
import java.time.Clock

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.SysEnv
import nextflow.exception.AbortSignalException
import org.fusesource.jansi.internal.CLibrary

/**
 * Tells the terminal how to resume this run with the Program Resume Protocol
 * (OSC 7577, revision 0.5), so that a terminal or multiplexer that restarts
 * can run {@code nextflow run ... -resume <sessionId>} again in the launch
 * directory and pick the pipeline up from the task cache.
 *
 * <p>One report is sent when the run is created, because the command never
 * changes during a run. When the run ends because the user or the pipeline
 * ended it, a {@code clear} report removes the record. When it ends because
 * the terminal went away (SIGHUP) or the system asked it to stop (SIGTERM),
 * nothing is sent: the spec forbids {@code clear} then, because that is the
 * moment the record is needed.
 *
 * <p>Off unless {@code NXF_OSC7577} is true, standard output is a terminal,
 * and the terminal supplied {@code OSC7577_TOKEN}.
 */
@Slf4j
@CompileStatic
class ResumeObserver implements TraceObserverV2 {

    static final String ENV = 'NXF_OSC7577'

    private final String token

    private final OutputStream out

    private final long pid

    private final Clock clock

    private final String cwd

    private Session session

    private boolean reported

    ResumeObserver(String token, OutputStream out, long pid, Clock clock, String cwd) {
        this.token = token
        this.out = out
        this.pid = pid
        this.clock = clock
        this.cwd = cwd
    }

    /**
     * @return an observer that writes to the terminal, or {@code null} when reporting is off
     */
    static ResumeObserver create() {
        if( !enabled(stdoutIsTty()) )
            return null
        return new ResumeObserver(
                SysEnv.get(Osc7577.TOKEN_ENV),
                new FileOutputStream(FileDescriptor.out),
                ProcessHandle.current().pid(),
                Clock.systemUTC(),
                Paths.get('').toAbsolutePath().toString() )
    }

    /**
     * Opt in with {@code NXF_OSC7577=true}. Without a token from the terminal there
     * is no way to build a report the terminal accepts, so that is also off.
     */
    static boolean enabled(boolean stdoutTty) {
        return SysEnv.getBool(ENV, false) && stdoutTty && SysEnv.get(Osc7577.TOKEN_ENV)
    }

    @Override
    void onFlowCreate(Session session) {
        this.session = session
        if( !session.resumeCommand ) {
            log.debug "OSC 7577: this run has no resume command"
            return
        }
        final text = Osc7577.report(token, pid, now(), session.resumeCommand, cwd)
        if( !text ) {
            log.debug "OSC 7577: the resume command cannot be reported -- ${session.resumeCommand}"
            return
        }
        reported = write(text)
    }

    @Override
    void onFlowComplete() {
        if( !reported || terminalEnded(session) )
            return
        final text = Osc7577.clear(token, pid, now())
        if( text )
            write(text)
    }

    /**
     * SIGHUP means the terminal went away and SIGTERM that the system is stopping
     * the run. A resume record is needed in both cases.
     */
    static boolean terminalEnded(Session session) {
        final error = session?.error
        return error instanceof AbortSignalException && error.message in ['SIGHUP', 'SIGTERM']
    }

    private long now() {
        return clock.millis().intdiv(1000)
    }

    private boolean write(String text) {
        try {
            out.write(text.getBytes('UTF-8'))
            out.flush()
            return true
        }
        catch( IOException e ) {
            log.debug "OSC 7577: cannot write to the terminal -- ${e.message}"
            return false
        }
    }

    /**
     * True when file descriptor 1 is a terminal. {@link System#console()} also needs a
     * console on standard input, and on recent JDKs it can be set without a tty, so
     * it is not used. Where the native check is missing, the answer is no.
     */
    static boolean stdoutIsTty() {
        try {
            return CLibrary.HAVE_ISATTY && CLibrary.isatty(CLibrary.STDOUT_FILENO) != 0
        }
        catch( Throwable ignored ) {
            return false
        }
    }

}
