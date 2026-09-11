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

package io.seqera.tower.plugin


import groovy.util.logging.Slf4j
import io.seqera.util.retry.Retryable
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.script.dsl.Description
import nextflow.util.Duration
import nextflow.util.RetryConfig

/**
 * Configuration class for Tower retry policy settings.
 *
 * This class defines the retry behavior for Tower operations including HTTP requests
 * and other potentially failing operations. It implements exponential backoff with
 * jitter to handle transient failures gracefully.
 *
 * The retry policy supports:
 * - Configurable initial delay before the first retry attempt
 * - Maximum delay cap to prevent excessively long wait times
 * - Limited number of retry attempts to avoid infinite loops
 * - Jitter randomization to prevent thundering herd problems
 * - Exponential backoff multiplier for progressive delay increases
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class TowerRetryPolicy implements Retryable.Config, ConfigScope {

    /**
     * Default max number of attempts. Combined with the default {@code delay} (350ms),
     * {@code multiplier} (2.0) and {@code maxDelay} (90s), the 9 exponential backoff gaps
     * span a retry window of about 3 minutes, so that transient gateway errors (e.g. a
     * `502 Bad Gateway`) can be ridden out before aborting the run.
     */
    static final int DEFAULT_MAX_ATTEMPTS = 10

    @ConfigOption
    @Description("""
        Initial delay before retrying a failed Tower operation (default: `350ms`).
    """)
    Duration delay

    @ConfigOption
    @Description("""
        Maximum delay between retry attempts for Tower operations (default: `90s`).
    """)
    Duration maxDelay

    @ConfigOption
    @Description("""
        Maximum number of attempts for Tower operations, including the initial one (default: `10`). Use `-1` for no limit.
    """)
    int maxAttempts

    @ConfigOption
    @Description("""
        Random jitter factor applied to retry delays to avoid thundering herd issues (default: `0.25`).
    """)
    double jitter

    @ConfigOption
    @Description("""
        Multiplier factor for exponential backoff between retry attempts (default: `2.0`).
    """)
    double multiplier

    TowerRetryPolicy(Map opts, Map legacy=Map.of()) {
        // note: use explicit null checks rather than the elvis operator, because 0 and 0.0
        // are falsy in Groovy -- `jitter = 0` or `maxAttempts = 0` are meaningful settings
        // and must not be silently replaced by the defaults
        this.delay = firstOf(opts.delay, legacy.backOffDelay, RetryConfig.DEFAULT_DELAY) as Duration
        this.maxDelay = firstOf(opts.maxDelay, RetryConfig.DEFAULT_MAX_DELAY) as Duration
        this.maxAttempts = firstOf(opts.maxAttempts, legacy.maxRetries, DEFAULT_MAX_ATTEMPTS) as Integer
        this.jitter = firstOf(opts.jitter, RetryConfig.DEFAULT_JITTER) as Double
        this.multiplier = firstOf(opts.multiplier, legacy.backOffBase, RetryConfig.DEFAULT_MULTIPLIER) as Double
        // `maxAttempts` counts the initial attempt, and -1 is the retry library's value for
        // "retry indefinitely". Anything else below 1 would mean "never run the operation at
        // all" and is rejected by the library, so warn and fall back to a single attempt --
        // which is what `0` was almost certainly meant to express -- rather than abort a run
        // over a setting that used to be silently ignored.
        if( maxAttempts < 1 && maxAttempts != -1 ) {
            log.warn "Invalid value for config option 'tower.retryPolicy.maxAttempts' -- offending value: $maxAttempts -- using 1 (no retries) instead"
            this.maxAttempts = 1
        }
    }

    /**
     * Return the first non-null value, so that a legitimate falsy setting such as
     * `0` or `0.0` is honoured instead of falling through to the default.
     */
    private static Object firstOf(Object... values) {
        for( Object it : values )
            if( it != null )
                return it
        return null
    }
}
