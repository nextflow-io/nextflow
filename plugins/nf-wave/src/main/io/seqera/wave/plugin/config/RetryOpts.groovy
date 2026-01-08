/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package io.seqera.wave.plugin.config

import groovy.transform.CompileStatic
import groovy.transform.ToString
import io.seqera.util.retry.Retryable
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.script.dsl.Description
import nextflow.util.Duration

/**
 * Model retry options for Wave http requests
 */
@ToString(includeNames = true, includePackage = false)
@CompileStatic
class RetryOpts implements ConfigScope, Retryable.Config {

    @ConfigOption
    @Description("""
        The initial delay when a failing HTTP request is retried (default: `450ms`).
    """)
    Duration delay = Duration.of('450ms')

    @ConfigOption
    @Description("""
        The max delay when a failing HTTP request is retried (default: `90s`).
    """)
    Duration maxDelay = Duration.of('90s')

    @ConfigOption
    @Description("""
        The max number of attempts a failing HTTP request is retried (default: `5`).
    """)
    int maxAttempts = 5

    @ConfigOption
    @Description("""
        The jitter factor used to randomly vary retry delays (default: `0.25`).
    """)
    double jitter = 0.25

    @ConfigOption
    @Description("""
        The multiplier used for exponential backoff delay calculations (default: `2.0`)
    """)
    double multiplier = 2;

    RetryOpts() {
        this(Collections.emptyMap())
    }

    RetryOpts(Map config) {
        if( config.delay )
            delay = config.delay as Duration
        if( config.maxDelay )
            maxDelay = config.maxDelay as Duration
        if( config.maxAttempts )
            maxAttempts = config.maxAttempts as int
        if( config.jitter )
            jitter = config.jitter as double
        if( config.multiplier )
            multiplier = config.multiplier as double
    }

    // Methods required by Retryable.Config interface
    @Override
    java.time.Duration getDelayAsDuration() {
        return java.time.Duration.ofMillis(delay.toMillis())
    }

    @Override
    java.time.Duration getMaxDelayAsDuration() {
        return java.time.Duration.ofMillis(maxDelay.toMillis())
    }
}
