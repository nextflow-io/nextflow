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

import nextflow.util.Duration
import nextflow.util.RetryConfig
import spock.lang.Specification

/**
 * Unit tests for TowerRetryPolicy
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TowerRetryPolicyTest extends Specification {

    def 'should validate default values of tower retry policy'() {
        when:
        def policy = new TowerRetryPolicy([:])
        
        then:
        policy.delay == RetryConfig.DEFAULT_DELAY
        policy.maxDelay == RetryConfig.DEFAULT_MAX_DELAY
        policy.maxAttempts == RetryConfig.DEFAULT_MAX_ATTEMPTS
        policy.jitter == RetryConfig.DEFAULT_JITTER
        policy.multiplier == RetryConfig.DEFAULT_MULTIPLIER
    }
    
    def 'should use provided values when specified'() {
        when:
        def customOptions = [
                delay: '1s' as nextflow.util.Duration,
                maxDelay: '60s' as nextflow.util.Duration,
                maxAttemps: 3,
                jitter: 0.5,
                multiplier: 1.5
        ]
        def policy = new TowerRetryPolicy(customOptions)
        
        then:
        policy.delay == customOptions.delay
        policy.maxDelay == customOptions.maxDelay
        policy.maxAttempts == 3
        policy.jitter == 0.5d
        policy.multiplier == 1.5d
    }

    def 'should use provided values when specified'() {
        when:
        def policy = new TowerRetryPolicy([:], [backOffDelay: 500, maxRetries: 100, backOffBase: 5])

        then:
        policy.delay == Duration.of('500ms')
        policy.maxAttempts == 100
        policy.multiplier == 5
        and:
        policy.maxDelay == RetryConfig.DEFAULT_MAX_DELAY
        policy.jitter == RetryConfig.DEFAULT_JITTER
    }
}
