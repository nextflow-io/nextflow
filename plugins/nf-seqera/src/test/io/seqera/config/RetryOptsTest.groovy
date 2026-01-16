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

import nextflow.util.Duration
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class RetryOptsTest extends Specification {

    def 'should create retry config' () {

        expect:
        new RetryOpts().delay0 == Duration.of('450ms')
        new RetryOpts().maxDelay0 == Duration.of('90s')
        new RetryOpts().maxAttempts0 == 10
        new RetryOpts().jitter0 == 0.25d
        new RetryOpts().multiplier0 == 2.0d

        and:
        new RetryOpts([maxAttempts: 20]).maxAttempts0 == 20
        new RetryOpts([delay: '1s']).delay0 == Duration.of('1s')
        new RetryOpts([maxDelay: '1m']).maxDelay0 == Duration.of('1m')
        new RetryOpts([jitter: '0.5']).jitter0 == 0.5d
        new RetryOpts([multiplier: '3.0']).multiplier0 == 3.0d

    }

    def 'should implement Retryable.Config interface' () {
        when:
        def opts = new RetryOpts([delay: '1s', maxDelay: '2m', maxAttempts: 5, jitter: '0.3', multiplier: '1.5'])

        then:
        opts.getDelay() == java.time.Duration.ofSeconds(1)
        opts.getMaxDelay() == java.time.Duration.ofMinutes(2)
        opts.getMaxAttempts() == 5
        opts.getJitter() == 0.3d
        opts.getMultiplier() == 1.5d
    }

}
