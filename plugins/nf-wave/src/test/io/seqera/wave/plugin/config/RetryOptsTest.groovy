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


import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class RetryOptsTest extends Specification {

    def 'should create retry config' () {
        expect:
        new RetryOpts().delay == java.time.Duration.ofMillis(450)
        new RetryOpts().maxDelay == java.time.Duration.ofSeconds(90)
        new RetryOpts().maxAttempts == 5
        new RetryOpts().jitter == 0.25d

        and:
        new RetryOpts([maxAttempts: 20]).maxAttempts == 20
        new RetryOpts([delay: '1s']).delay == java.time.Duration.ofSeconds(1)
        new RetryOpts([maxDelay: '1m']).maxDelay == java.time.Duration.ofMinutes(1)
        new RetryOpts([jitter: '0.5']).jitter == 0.5d
    }

}
