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

package nextflow.util

import nextflow.SysEnv
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class RetryConfigTest extends Specification {

    def 'should create retry config' () {
        expect:
        new RetryConfig().delay == java.time.Duration.ofMillis(350)
        new RetryConfig().maxDelay == java.time.Duration.ofSeconds(90)
        new RetryConfig().maxAttempts == 5
        new RetryConfig().jitter == 0.25d
        new RetryConfig().multiplier == 2d

        and:
        new RetryConfig([maxAttempts: 20]).maxAttempts == 20
        new RetryConfig([delay: '1s']).delay == java.time.Duration.ofSeconds(1)
        new RetryConfig([maxDelay: '1m']).maxDelay == java.time.Duration.ofMinutes(1)
        new RetryConfig([jitter: '0.5']).jitter == 0.5d
        new RetryConfig([multiplier: '5']).multiplier == 5d
    }

    def 'should get the setting from the system env' () {
        given:
        SysEnv.push([
            NXF_RETRY_POLICY_DELAY: '10s',
            NXF_RETRY_POLICY_MAX_DELAY: '100s',
            NXF_RETRY_POLICY_MAX_ATTEMPTS: '1000',
            NXF_RETRY_POLICY_JITTER: '10000',
            NXF_RETRY_POLICY_MULTIPLIER: '90'
        ])

        expect:
        new RetryConfig().getDelay() == java.time.Duration.ofSeconds(10)
        new RetryConfig().getMaxDelay() == java.time.Duration.ofSeconds(100)
        new RetryConfig().getMaxAttempts() == 1000
        new RetryConfig().getJitter() == 10_000
        new RetryConfig().getMultiplier() == 90d

        cleanup:
        SysEnv.pop()
    }

    @Unroll
    def 'should get config from map' () {
        given:
        def NAME = 'foo'
        def PREFIX = 'P_'
        and:
        SysEnv.push(ENV)

        expect:
        RetryConfig.valueOf(CONFIG, NAME, PREFIX, DEF_VAL, DEF_TYPE) == EXPECTED

        cleanup:
        SysEnv.pop()

        where:
        CONFIG          | ENV           | DEF_VAL       | DEF_TYPE      | EXPECTED
        null            | [:]           | null          | String        | null
        [:]             | [:]           | null          | String        | null
        [:]             | [:]           | 'one'         | String        | 'one'
        [foo:'two']     | [:]           | 'one'         | String        | 'two'
        [foo:'']        | [:]           | 'one'         | String        | ''
        [foo:'two']     | [P_FOO:'bar'] | 'one'         | String        | 'two'
        [:]             | [P_FOO:'bar'] | 'one'         | String        | 'bar'

        and:
        null            | [:]           | null          | Integer       | null
        [:]             | [:]           | null          | Integer       | null
        [:]             | [:]           | 1             | Integer       | 1
        [foo:2]         | [:]           | 1             | Integer       | 2
        [foo:'2']       | [:]           | 1             | Integer       | 2
        [foo:'2']       | [P_FOO:'3']   | 1             | Integer       | 2
        [:]             | [P_FOO:'3']   | 1             | Integer       | 3

        and:
        null            | [:]           | null          | Boolean       | null
        [:]             | [:]           | true          | Boolean       | true
        [foo:false]     | [:]           | true          | Boolean       | false
        [foo:'false']   | [:]           | true          | Boolean       | false
        [foo:true]      | [:]           | false         | Boolean       | true
        [foo:'true']    | [:]           | false         | Boolean       | true
        [foo:'true']    | [P_FOO:'false']| null         | Boolean       | true
        [:]             | [P_FOO:'false']| null         | Boolean       | false
        [:]             | [P_FOO:'true'] | null         | Boolean       | true

        and:
        [:]             | [:]           | Duration.of('1s') | Duration  | Duration.of('1s')
        [foo:'10ms']    | [:]           | null              | Duration  | Duration.of('10ms')
        [:]             | [P_FOO:'1s']  | null              | Duration  | Duration.of('1s')
    }

    def 'should map camelCase to snake uppercase' () {
        given:
        SysEnv.push(ENV)

        expect:
        RetryConfig.valueOf([:], NAME, PREFIX, null, String) == EXPECTED

        cleanup:
        SysEnv.pop()

        where:
        EXPECTED    | PREFIX    | NAME              | ENV
        null        | 'foo'     | 'bar'             | [:]
        'one'       | 'foo'     | 'bar'             | [FOO_BAR: 'one']
        'one'       | 'foo_'    | 'bar'             | [FOO_BAR: 'one']
        'one'       | 'foo_'    | 'thisAndThat'     | [FOO_THIS_AND_THAT: 'one']
    }
}
