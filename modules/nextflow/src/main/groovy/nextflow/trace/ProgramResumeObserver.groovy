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
import java.util.function.LongSupplier

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.SysEnv
import nextflow.exception.AbortSignalException

/**
 * Tells the terminal how to resume this run, using the Program Resume
 * Protocol (OSC 7577).
 *
 * <p>One report when the run starts, naming {@code nextflow run ... -resume <session id>}
 * and the launch directory. One {@code clear} when the user is finished with
 * the run: it completed, failed, or was interrupted with Ctrl-C. No
 * {@code clear} when the run ends because the terminal went away (SIGHUP or
 * SIGTERM), which is when the record is needed.
 *
 * <p>Nothing is written unless the run is opted in (see {@link ProgramResume}),
 * the terminal supplied {@code OSC7577_TOKEN}, and stdout is a terminal.
 */
@Slf4j
@CompileStatic
class ProgramResumeObserver implements TraceObserverV2 {

    private static final Set<String> TERMINAL_ENDED = ['TERM', 'HUP'] as Set

    private final Appendable out

    private final boolean tty

    private final long pid

    private final LongSupplier epochSeconds

    private Session session

    private Osc7577Emitter emitter

    ProgramResumeObserver() {
        this(System.out, Osc7501Emitter.stdoutIsTty(), ProcessHandle.current().pid(), { System.currentTimeMillis().intdiv(1000) as long })
    }

    ProgramResumeObserver(Appendable out, boolean tty, long pid, LongSupplier epochSeconds) {
        this.out = out
        this.tty = tty
        this.pid = pid
        this.epochSeconds = epochSeconds
    }

    @Override
    void onFlowCreate(Session session) {
        this.session = session
        final skip = skipReason(session)
        if( skip ) {
            log.debug "OSC 7577 resume not reported: ${skip}"
            return
        }
        final cwd = Paths.get('').toAbsolutePath()
        final argv = ProgramResume.command(ProgramResume.launcher(SysEnv.get('NXF_CLI'), cwd), session.cliArgs, session.uniqueId.toString())
        if( argv == null ) {
            log.debug 'OSC 7577 resume not reported: not a run command'
            return
        }
        try {
            final candidate = new Osc7577Emitter(out, ProgramResume.token(), pid, epochSeconds)
            candidate.report(argv, cwd.toString())
            emitter = candidate
            log.debug "OSC 7577 resume reported: ${argv.join(' ')}"
        }
        catch( IllegalArgumentException e ) {
            log.debug "OSC 7577 resume not reported: ${e.message}"
        }
    }

    @Override
    void onFlowComplete() {
        if( !emitter || terminalEnded(session) )
            return
        emitter.clear()
    }

    private String skipReason(Session session) {
        if( !tty )
            return 'stdout is not a terminal'
        if( !Osc7577Emitter.isValidToken(ProgramResume.token()) )
            return "${ProgramResume.TOKEN_ENV} is not set to a valid token"
        if( session.preview )
            return 'preview run'
        return null
    }

    /**
     * The terminal ending sends SIGHUP, and a system shutdown sends SIGTERM. The
     * protocol forbids clearing the record then. Ctrl-C is the user ending the run.
     */
    static boolean terminalEnded(Session session) {
        final error = session?.error
        return error instanceof AbortSignalException && ((AbortSignalException) error).signalName in TERMINAL_ENDED
    }

}
