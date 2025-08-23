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

package io.seqera.tower.plugin


import io.seqera.util.retry.Retryable
import nextflow.config.schema.ConfigOption
import nextflow.config.schema.ConfigScope
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
class TowerRetryPolicy implements Retryable.Config, ConfigScope {

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
    Maximum number of retry attempts for Tower operations (default: `5`).
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
        this.delay = opts.delay as Duration ?: legacy.backOffDelay as Duration ?: RetryConfig.DEFAULT_DELAY
        this.maxDelay = opts.maxDelay as Duration ?: RetryConfig.DEFAULT_MAX_DELAY
        this.maxAttempts = opts.maxAttemps as Integer ?: legacy.maxRetries as Integer ?: RetryConfig.DEFAULT_MAX_ATTEMPTS
        this.jitter = opts.jitter as Double ?: RetryConfig.DEFAULT_JITTER
        this.multiplier = opts.multiplier as Double ?: legacy.backOffBase as Double ?: RetryConfig.DEFAULT_MULTIPLIER
    }
}
