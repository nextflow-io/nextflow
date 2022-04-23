/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.cloud.azure.config

import nextflow.util.Duration
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzRetryConfigTest extends Specification {

    def 'should create retry config' () {

        expect:
        new AzRetryConfig().delay == Duration.of('50ms')
        new AzRetryConfig().maxDelay == Duration.of('30s')
        new AzRetryConfig().maxAttempts == 5
        new AzRetryConfig().jitter == 0.25d

        and:
        new AzRetryConfig([maxAttempts: 10]).maxAttempts == 10
        new AzRetryConfig([delay: '1s']).delay == Duration.of('1s')
        new AzRetryConfig([maxDelay: '1m']).maxDelay == Duration.of('1m')
        new AzRetryConfig([jitter: '0.5']).jitter == 0.5d

    }
}
