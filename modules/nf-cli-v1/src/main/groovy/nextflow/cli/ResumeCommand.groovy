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

import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.CompileStatic

/**
 * Builds the command line that resumes a `nextflow run` as the same session.
 *
 * <p>The command is the original one with its run options kept, its script
 * arguments untouched, any {@code -resume} replaced by {@code -resume <sessionId>}
 * and any {@code -name} dropped (a resumed run gets a new name, and reusing the
 * old one is an error). Each element stays a separate argument, so nothing is
 * ever parsed by a shell.
 */
@CompileStatic
class ResumeCommand {

    static final String DEFAULT_LAUNCHER = 'nextflow'

    /**
     * @param launcher the program to execute, see {@link #launcherPath}
     * @param args the arguments as {@link Launcher#getNormalizedArgs}, where `-resume` is always followed by its value
     * @param sessionId the session to resume
     * @return the command, or {@code null} when it cannot be resumed: it is not a `run`, or the script comes from standard input
     */
    static List<String> of(String launcher, List<String> args, UUID sessionId) {
        final result = new ArrayList<String>(args.size() + 3)
        for( int i = 0; i < args.size(); i++ ) {
            final arg = args[i]
            if( arg == '-resume' || arg == '-name' ) {
                // skip the value too
                i++
                continue
            }
            result.add(arg)
        }
        final run = result.indexOf('run')
        if( run == -1 )
            return null
        if( run + 1 < result.size() && result[run + 1] == '-' )
            return null
        result.add(run + 1, '-resume')
        result.add(run + 2, sessionId.toString())
        result.add(0, launcher)
        return result
    }

    /**
     * The launcher script exports {@code NXF_CLI="$0 $@"}, so what precedes the
     * arguments is the path the user ran. A relative path is made absolute against
     * the launch directory, because the terminal restores in another one and the
     * spec forbids relative paths. A bare name stays bare and is looked up on
     * PATH when the command is resumed.
     *
     * @param nxfCli the value of {@code NXF_CLI}
     * @param args the arguments exactly as the user typed them
     * @return the program to execute, {@code nextflow} when it cannot be told
     */
    static String launcherPath(String nxfCli, String[] args, Path cwd = Paths.get('').toAbsolutePath()) {
        if( !nxfCli )
            return DEFAULT_LAUNCHER
        final suffix = ' ' + args.join(' ')
        if( !nxfCli.endsWith(suffix) )
            return DEFAULT_LAUNCHER
        final path = nxfCli.substring(0, nxfCli.length() - suffix.length())
        if( !path )
            return DEFAULT_LAUNCHER
        return path.contains('/') ? cwd.resolve(path).normalize().toString() : path
    }

}
