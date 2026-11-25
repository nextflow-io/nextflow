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

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.SysEnv

/**
 * Decides whether a run takes part in the Program Resume Protocol (OSC 7577)
 * and what resume command it reports.
 *
 * <p>A run reports only when it is opted in with {@code NXF_OSC7577} or the
 * {@code osc7577} config option. The default is off, because a terminal that
 * accepts the report will start the command again without the user typing it.
 *
 * <p>The reported command is the launch command with {@code -resume <session id>}
 * in place of any {@code -resume} or {@code -name} the user gave, so the
 * terminal re-launches the same run on top of the task cache.
 */
@CompileStatic
class ProgramResume {

    static final String ENV = 'NXF_OSC7577'

    static final String CONFIG = 'osc7577'

    /** Set by the terminal in the first process of each pseudo-terminal; its presence is how support is detected */
    static final String TOKEN_ENV = 'OSC7577_TOKEN'

    private static final String RUN = 'run'

    private static final String LAUNCHER = 'nextflow'

    /**
     * Opt-in only. {@code NXF_OSC7577}, when set to a non-blank value, wins over
     * the config option, so an environment variable can switch a configured run off.
     */
    static boolean optedIn(Session session) {
        final raw = SysEnv.get(ENV)
        if( raw?.trim() )
            return truthy(raw)
        final config = session?.config
        if( config == null )
            return false
        final value = config.get(CONFIG)
        return value != null && truthy(value.toString())
    }

    static String token() {
        return SysEnv.get(TOKEN_ENV)
    }

    /**
     * Build the command that brings a run back.
     *
     * @param launcher the program to execute, a bare name or an absolute path
     * @param args the normalized launcher arguments, starting with the top-level options and the {@code run} command
     * @param sessionId the session unique id, which names the task cache to resume
     * @return the argv, or {@code null} when {@code args} is not a {@code run} command
     */
    static List<String> command(String launcher, List<String> args, String sessionId) {
        final run = args ? args.indexOf(RUN) : -1
        if( run < 0 )
            return null
        final result = new ArrayList<String>(args.size() + 3)
        result.add(launcher)
        result.addAll(args.subList(0, run + 1))
        result.add('-resume')
        result.add(sessionId)
        // a rerun cannot reuse the run name, and a second -resume would take the place of ours
        int i = run + 1
        while( i < args.size() ) {
            final arg = args[i++]
            if( arg == '-name' || arg == '-resume' )
                i++
            else
                result.add(arg)
        }
        return result
    }

    /**
     * The launcher path as the user typed it, made absolute, so the resumed run
     * uses the same {@code nextflow} even if {@code PATH} differs after a restart.
     * The launcher script exports {@code NXF_CLI} as {@code "$0 $@"}. When
     * {@code $0} cannot be told apart from the arguments, fall back to the bare
     * name and let the terminal look it up on {@code PATH}.
     */
    static String launcher(String nxfCli, Path cwd) {
        if( !nxfCli )
            return LAUNCHER
        final space = nxfCli.indexOf(' ')
        final program = space < 0 ? nxfCli : nxfCli.substring(0, space)
        if( program == LAUNCHER || !program.endsWith('/' + LAUNCHER) )
            return LAUNCHER
        return cwd.resolve(program).normalize().toString()
    }

    /**
     * A program that reports must not hand the token to what it starts, or a
     * task could write its own resume record over the run's. Only runs that
     * are opted in lose the variable.
     */
    static void withholdToken(ProcessBuilder builder, Session session) {
        if( optedIn(session) )
            builder.environment().remove(TOKEN_ENV)
    }

    private static boolean truthy(String text) {
        final value = text.trim()
        return Boolean.parseBoolean(value) || value == '1'
    }

}
