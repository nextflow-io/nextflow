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

package nextflow.util

import nextflow.SysEnv
import spock.lang.Specification
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class RetryConfigTest extends Specification {

    def 'should create retry config' () {

        expect:
        new RetryConfig().delay == Duration.of('350ms')
        new RetryConfig().maxDelay == Duration.of('90s')
        new RetryConfig().maxAttempts == 5
        new RetryConfig().jitter == 0.25d

        and:
        new RetryConfig([maxAttempts: 20]).maxAttempts == 20
        new RetryConfig([delay: '1s']).delay == Duration.of('1s')
        new RetryConfig([maxDelay: '1m']).maxDelay == Duration.of('1m')
        new RetryConfig([jitter: '0.5']).jitter == 0.5d

    }

    def 'should get the setting from the system env' () {
        given:
        SysEnv.push([
            NXF_RETRY_POLICY_DELAY: '10s',
            NXF_RETRY_POLICY_MAX_DELAY: '100s',
            NXF_RETRY_POLICY_MAX_ATTEMPTS: '1000',
            NXF_RETRY_POLICY_JITTER: '10000'
        ])

        expect:
        new RetryConfig().getDelay() == Duration.of('10s')
        new RetryConfig().getMaxDelay() == Duration.of('100s')
        new RetryConfig().getMaxAttempts() == 1000
        new RetryConfig().getJitter() == 10_000

        cleanup:
        SysEnv.pop()
    }

}
