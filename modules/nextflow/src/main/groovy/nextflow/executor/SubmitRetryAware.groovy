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
import nextflow.processor.TaskRun
import nextflow.util.Duration
/**
 * Generic retry-able submit implementation for executors.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
trait SubmitRetryAware {

    abstract Executor getExecutor()

    abstract TaskRun getTask()

    <T> T safeExecute(CheckedSupplier<T> action) {
        Failsafe.with(retryPolicy()).get(action)
    }

    private <T> RetryPolicy<T> retryPolicy() {

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
    private Predicate<? extends Throwable> retryCondition(String reasonPattern) {
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
