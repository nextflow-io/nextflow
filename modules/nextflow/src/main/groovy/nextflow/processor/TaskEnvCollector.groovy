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

package nextflow.processor

import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.time.temporal.ChronoUnit
import java.util.regex.Matcher

import dev.failsafe.Failsafe
import dev.failsafe.FailsafeException
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessEvalException
import nextflow.util.RetryConfig
/**
 * Implements the collection of environment variables
 * from the environment of a task execution.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskEnvCollector {

    private Path workDir

    private Map<String,String> evalCmds

    private RetryConfig retryConfig

    TaskEnvCollector(Path workDir, Map<String,String> evalCmds) {
        this(workDir, evalCmds, RetryConfig.config())
    }

    protected TaskEnvCollector(Path workDir, Map<String,String> evalCmds, RetryConfig retryConfig) {
        this.workDir = workDir
        this.evalCmds = evalCmds
        this.retryConfig = retryConfig
    }

    /**
     * Load the values for `env` and `eval` outputs from the `.command.env` file.
     *
     * The file is written by the task wrapper script and may not be immediately
     * visible once the task has completed (e.g. due to propagation delays on
     * shared/eventually-consistent filesystems), so missing or empty content is
     * retried for a short, bounded window before giving up.
     */
    Map collect() {
        final env = readEnvFile()
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

    /**
     * Read the `.command.env` file, retrying when it's missing or empty.
     */
    private String readEnvFile() {
        final envFile = workDir.resolve(TaskRun.CMD_ENV)
        final listener = new EventListener<ExecutionAttemptedEvent>() {
            @Override
            void accept(ExecutionAttemptedEvent event) throws Throwable {
                log.debug "Failed to read task env file: ${envFile.toUriString()} -- attempt: ${event.attemptCount}; reason: ${event.lastException.message}"
            }
        }
        final retryPolicy = RetryPolicy.builder()
            .handle(NoSuchFileException, EmptyEnvFileException)
            .withBackoff(retryConfig.delay.toMillis(), retryConfig.maxDelay.toMillis(), ChronoUnit.MILLIS)
            .withMaxAttempts(retryConfig.maxAttempts)
            .withJitter(retryConfig.jitter)
            .onRetry(listener)
            .build()
        try {
            return Failsafe
                .with(retryPolicy)
                .get({it -> readEnvFile0(envFile)})
        }
        catch( FailsafeException e ) {
            // unwrap to surface the original missing/empty-file error instead of
            // the generic failsafe wrapper
            throw e.cause ?: e
        }
    }

    private String readEnvFile0(Path envFile) {
        final text = envFile.text
        if( !text )
            throw new EmptyEnvFileException(envFile)
        return text
    }

    static private class EmptyEnvFileException extends RuntimeException {
        EmptyEnvFileException(Path envFile) {
            super("Task env file is empty: ${envFile.toUriString()}".toString())
        }
    }
}
