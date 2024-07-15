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

package nextflow.cloud.azure.batch

import java.time.Duration

import com.azure.storage.common.policy.RetryPolicyType
import nextflow.cloud.azure.config.AzRetryConfig
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzHelperTest extends Specification {

    @Unroll
    def 'should create retry options' () {
        given:
        def opts = AzHelper.requestRetryOptions0(CONFIG)
        expect:
        opts.@retryPolicyType == RetryPolicyType.EXPONENTIAL
        opts.maxRetryDelay.toMillis() == CONFIG.maxDelay.millis
        opts.retryDelay.toMillis() == CONFIG.delay.millis
        opts.maxTries == CONFIG.maxAttempts
        and:
        opts.secondaryHost == null
        opts.tryTimeoutDuration == Duration.ofSeconds(Integer.MAX_VALUE)

        where:
        _ | CONFIG
        _ | new AzRetryConfig()
        _ | new AzRetryConfig([delay: '10s', maxAttempts: 20, maxDelay: '200s'])
    }

}
