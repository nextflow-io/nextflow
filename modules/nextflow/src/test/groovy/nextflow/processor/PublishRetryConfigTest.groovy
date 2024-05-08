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

package nextflow.processor

import nextflow.util.Duration
import spock.lang.Specification

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class PublishRetryConfigTest extends Specification {

    def 'should create retry config' () {

        expect:
        new PublishRetryConfig().delay == Duration.of('350ms')
        new PublishRetryConfig().maxDelay == Duration.of('90s')
        new PublishRetryConfig().maxAttempts == 5
        new PublishRetryConfig().jitter == 0.25d

        and:
        new PublishRetryConfig([maxAttempts: 20]).maxAttempts == 20
        new PublishRetryConfig([delay: '1s']).delay == Duration.of('1s')
        new PublishRetryConfig([maxDelay: '1m']).maxDelay == Duration.of('1m')
        new PublishRetryConfig([jitter: '0.5']).jitter == 0.5d

    }

}
