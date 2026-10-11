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

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.Session
import nextflow.SysEnv

/**
 * Maps a Nextflow run onto OSC 7501 states and emits only when the state changes.
 *
 * <p>{@code working} while any task is pending, submitted, or running.
 * {@code idle} when the run is still going and nothing is in flight, and when
 * the user cancels or interrupts (the protocol asks for {@code idle} then).
 * {@code done} when the run completes successfully. {@code error} when it
 * stops because of a failure. {@code blocked} is not used: a run does not
 * pause for permission, a question, or a credential.
 */
@CompileStatic
class TerminalActivityTracker {

    static final String ENV = 'NXF_OSC7501'

    static final String CONFIG = 'osc7501'

    private final TerminalActivityEmitter emitter

    private String current

    private boolean settled

    private boolean suspended

    private TerminalActivityTracker(TerminalActivityEmitter emitter) {
        this.emitter = emitter
    }

    /**
     * @return a tracker that has already reported {@code idle}, or {@code null}
     *         when reporting is off
     */
    static TerminalActivityTracker start(Session session) {
        if( !enabled(session, Osc7501Emitter.stdoutIsTty()) )
            return null
        final tracker = new TerminalActivityTracker(new Osc7501Emitter())
        tracker.report('idle', 'Waiting for tasks')
        return tracker
    }

    @PackageScope
    static TerminalActivityTracker open(TerminalActivityEmitter emitter) {
        return new TerminalActivityTracker(emitter)
    }

    /**
     * Reporting is on only when ANSI logging is on, stdout is a terminal, and
     * neither {@code NXF_OSC7501} nor {@code osc7501} opts out. An unset
     * switch leaves the default (on) in place. A set switch that is not a
     * truthy boolean turns reporting off.
     */
    static boolean enabled(Session session, boolean stdoutTty) {
        if( !session?.ansiLog || !stdoutTty )
            return false
        return envAllows() && configAllows(session.config)
    }

    synchronized void onStats(int inFlight) {
        if( suspended || settled )
            return
        if( inFlight > 0 )
            report('working', 'Running tasks')
        else
            report('idle', 'Waiting for tasks')
    }

    /**
     * Ignore further queue updates. Called once the run is finishing so the
     * last progress paint does not report {@code idle} ahead of {@code done}
     * or {@code error}.
     */
    synchronized void beginShutdown() {
        suspended = true
    }

    synchronized void onFlowError(boolean interrupted) {
        if( settled )
            return
        settled = true
        if( interrupted )
            report('idle', 'Run cancelled')
        else
            report('error', 'Workflow failed')
    }

    synchronized void onRunComplete(boolean success, boolean interrupted) {
        if( settled )
            return
        settled = true
        if( interrupted )
            report('idle', 'Run cancelled')
        else if( success )
            report('done', 'Workflow completed')
        else
            report('error', 'Workflow failed')
    }

    private void report(String state, String message) {
        if( state == current )
            return
        current = state
        try {
            emitter.emit(state, message)
        }
        catch( IOException ignored ) {
            // a closed terminal must not fail the run, and this path renders
            // logs, so it cannot log the failure either
        }
    }

    private static boolean envAllows() {
        final raw = SysEnv.get(ENV)
        if( raw == null )
            return true
        final text = raw.trim()
        if( !text )
            return true
        return truthy(text)
    }

    private static boolean configAllows(Map config) {
        if( config == null || !config.containsKey(CONFIG) )
            return true
        final value = config.get(CONFIG)
        if( value == null )
            return true
        if( value instanceof Boolean )
            return (Boolean) value
        final text = value.toString().trim()
        if( !text )
            return true
        return truthy(text)
    }

    private static boolean truthy(String text) {
        return Boolean.parseBoolean(text) || text == '1'
    }

}
