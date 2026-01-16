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

package io.seqera.config

import java.time.temporal.TemporalAmount

import groovy.transform.CompileStatic
import groovy.transform.ToString
import io.seqera.util.retry.Retryable
import nextflow.util.Duration

/**
 * Model retry options for Seqera scheduler HTTP requests.
 * Implements {@link Retryable.Config} for integration with lib-retry.
 */
@ToString(includeNames = true, includePackage = false)
@CompileStatic
class RetryOpts implements Retryable.Config {
    Duration delay0 = Duration.of('450ms')
    Duration maxDelay0 = Duration.of('90s')
    int maxAttempts0 = 10
    double jitter0 = 0.25
    double multiplier0 = 2.0d

    RetryOpts() {
        this(Collections.emptyMap())
    }

    RetryOpts(Map config) {
        if( config.delay )
            delay0 = config.delay as Duration
        if( config.maxDelay )
            maxDelay0 = config.maxDelay as Duration
        if( config.maxAttempts )
            maxAttempts0 = config.maxAttempts as int
        if( config.jitter )
            jitter0 = config.jitter as double
        if( config.multiplier )
            multiplier0 = config.multiplier as double
    }

    @Override
    TemporalAmount getDelay() {
        return java.time.Duration.ofMillis(delay0.toMillis())
    }

    @Override
    TemporalAmount getMaxDelay() {
        return java.time.Duration.ofMillis(maxDelay0.toMillis())
    }

    @Override
    int getMaxAttempts() {
        return maxAttempts0
    }

    @Override
    double getJitter() {
        return jitter0
    }

    @Override
    double getMultiplier() {
        return multiplier0
    }
}
