/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.processor

import java.nio.file.Path
import java.util.regex.Matcher

import groovy.transform.CompileStatic
import nextflow.exception.ProcessEvalException
/**
 * Implements the collection of environment variables
 * from the environment of a task execution.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class TaskEnvCollector {

    private Path workDir

    private Map<String,String> evalCmds

    TaskEnvCollector(Path workDir, Map<String,String> evalCmds) {
        this.workDir = workDir
        this.evalCmds = evalCmds
    }

    /**
     * Load the values for `env` and `eval` outputs from the `.command.env` file.
     */
    Map collect() {
        final env = workDir.resolve(TaskRun.CMD_ENV).text
        final result = new HashMap<String,String>(50)
        Matcher matcher
        // `current` represents the current capturing env variable name
        String current = null
        for( String line : env.readLines() ) {
            // Opening condition:
            // line should match a KEY=VALUE syntax
            if( !current && (matcher = (line=~/([a-zA-Z_][a-zA-Z0-9_]*)=(.*)/)) ) {
                final key = matcher.group(1)
                final value = matcher.group(2)
                if (!key) continue
                result.put(key, value)
                current = key
            }
            // Closing condition:
            // line should match /KEY/ or /KEY/=exit_status
            else if( current && (matcher = (line=~/\/${current}\/(?:=exit:(\d+))?/)) ) {
                final status = matcher.group(1) as Integer ?: 0
                // when exit status is defined and it is a non-zero, it should be interpreted
                // as a failure of the execution of the output command; in this case the variable
                // holds the std error message
                if( evalCmds != null && status ) {
                    final cmd = evalCmds.get(current)
                    final out = result[current]
                    throw new ProcessEvalException("Unable to evaluate output", cmd, out, status)
                }
                // reset current key
                current = null
            }
            else if( current && line != null ) {
                result[current] += '\n' + line
            }
        }
        return result
    }
}
