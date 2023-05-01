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

package nextflow.executor

import java.time.temporal.ChronoUnit
import java.util.function.Predicate
import java.util.regex.Pattern

import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedSupplier
import groovy.transform.Memoized
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessNonZeroExitStatusException
import nextflow.fusion.FusionAwareTask
import nextflow.processor.TaskRun
import nextflow.util.CmdLineHelper
import nextflow.util.Duration
import org.slf4j.Logger
import org.slf4j.LoggerFactory
/**
 * Implementation of job submission for grid executors.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
trait SubmitJobAware extends FusionAwareTask {

    static private Logger log = LoggerFactory.getLogger(SubmitJobAware)

    abstract AbstractGridExecutor getExecutor()

    abstract TaskRun getTask()

    ProcessBuilder createProcessBuilder(boolean pipeLauncherScript) {

        // -- log the submit command line
        final cli = executor.getSubmitCommandLine(task, task.workDir.resolve(TaskRun.CMD_RUN), pipeLauncherScript)
        log.trace "start process ${task.name} > cli: ${cli}"

        // -- create the submit process
        ProcessBuilder builder = new ProcessBuilder()
            .command( cli as String[] )
            .redirectErrorStream(true)

        if( !pipeLauncherScript )
            builder .directory(task.workDir.toFile())

        return builder
    }

    String processStart(ProcessBuilder builder, String pipeScript) {
        final process = builder.start()

        try {
            // -- forward the job launcher script to the command stdin if required
            if( pipeScript ) {
                log.trace "[${executor.name.toUpperCase()}] Submit STDIN command ${task.name} >\n${pipeScript.indent()}"
                process.out << pipeScript
                process.out.close()
            }

            // -- wait the the process completes
            final result = process.text
            final exitStatus = process.waitFor()
            final cmd = launchCmd0(builder,pipeScript)

            if( exitStatus ) {
                throw new ProcessNonZeroExitStatusException("Failed to submit process to grid scheduler for execution", result, exitStatus, cmd)
            }

            // -- return the process stdout
            return result
        }
        finally {
            // make sure to release all resources
            process.in.closeQuietly()
            process.out.closeQuietly()
            process.err.closeQuietly()
            process.destroy()
        }
    }

    private String launchCmd0(ProcessBuilder builder, String pipeScript) {
        def result = CmdLineHelper.toLine(builder.command())
        if( pipeScript ) {
            result = "cat << 'LAUNCH_COMMAND_EOF' | ${result}\n"
            result += pipeScript.trim() + '\n'
            result += 'LAUNCH_COMMAND_EOF\n'
        }
        return result
    }

    <T> T safeExecute(CheckedSupplier<T> action) {
        Failsafe.with(retryPolicy()).get(action)
    }

    <T> RetryPolicy<T> retryPolicy() {

        final delay = executor.session.getConfigAttribute('executor.retry.delay', '500ms') as Duration
        final maxDelay = executor.session.getConfigAttribute('executor.retry.maxDelay', '30s') as Duration
        final jitter = executor.session.getConfigAttribute('executor.retry.jitter', '0.25') as double
        final maxAttempts = executor.session.getConfigAttribute('executor.retry.maxAttempts', '3') as int
        final reason = executor.session.getConfigAttribute('executor.submit.retry.reason', 'Socket timed out') as String

        RetryPolicy.<T>builder()
                .handleIf(retryCondition(reason))
                .withBackoff(delay.toMillis(), maxDelay.toMillis(), ChronoUnit.MILLIS)
                .withMaxAttempts(maxAttempts)
                .withJitter(jitter)
                .onFailedAttempt(new RetryListener(task: task))
                .build()
    }

    static private class RetryListener extends EventListener<ExecutionAttemptedEvent> {
        TaskRun task

        @Override
        void accept(ExecutionAttemptedEvent event) throws Throwable {
            final failure = event.getLastFailure()
            if( failure instanceof ProcessNonZeroExitStatusException ) {
                final msg = """\
                    Failed to submit process '${task.name}'
                        - attempt : ${event.attemptCount}
                        - command : ${failure.command}
                        - reason  : ${failure.reason}
                    """.stripIndent(true)
                log.warn msg
            } else {
                log.debug("Unexpected retry failure: ${failure?.message}", failure)
            }
        }
    }

    @Memoized
    Predicate<? extends Throwable> retryCondition(String reasonPattern) {
        new RetryPredicate(pattern: Pattern.compile(reasonPattern))
    }

    static private class RetryPredicate extends Predicate<Throwable> {
        Pattern pattern

        @Override
        boolean test(Throwable failure) {
            if( failure instanceof ProcessNonZeroExitStatusException ) {
                final reason = failure.reason
                return reason ? pattern.matcher(reason).find() : false
            }
            return false
        }
    }

}
