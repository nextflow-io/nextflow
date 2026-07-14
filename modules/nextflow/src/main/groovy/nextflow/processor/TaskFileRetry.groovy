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

import java.time.temporal.ChronoUnit

import dev.failsafe.Failsafe
import dev.failsafe.FailsafeException
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.util.RetryConfig
/**
 * Generic retry wrapper for reading task-produced files (e.g. {@code .exitcode},
 * {@code .command.env}, {@code .command.out}) that may not be immediately visible
 * once a task has completed, due to propagation delays on shared or
 * eventually-consistent filesystems.
 *
 * Callers decide which exceptions should trigger a retry, since "missing" and
 * "empty" don't mean the same thing for every file (e.g. an empty task env file
 * always indicates an incomplete write, while empty stdout can be a legitimate
 * result).
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskFileRetry {

    static <T> T withRetry(RetryConfig retryConfig, List<Class<? extends Throwable>> retryOn, String label, Closure<T> body) {
        final listener = new EventListener<ExecutionAttemptedEvent>() {
            @Override
            void accept(ExecutionAttemptedEvent event) throws Throwable {
                log.debug "Failed to read ${label} -- attempt: ${event.attemptCount}; reason: ${event.lastException.message}"
            }
        }
        final retryPolicy = RetryPolicy.builder()
            .handle(retryOn)
            .withBackoff(retryConfig.delay.toMillis(), retryConfig.maxDelay.toMillis(), ChronoUnit.MILLIS)
            .withMaxAttempts(retryConfig.maxAttempts)
            .withJitter(retryConfig.jitter)
            .onRetry(listener)
            .build()
        try {
            return Failsafe
                .with(retryPolicy)
                .get({ it -> body.call() })
        }
        catch( FailsafeException e ) {
            // unwrap to surface the original missing/empty-file error instead of
            // the generic failsafe wrapper
            throw e.cause ?: e
        }
    }
}
