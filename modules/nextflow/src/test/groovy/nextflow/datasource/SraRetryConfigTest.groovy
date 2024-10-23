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

package nextflow.datasource

import nextflow.util.Duration
import spock.lang.Specification

/**
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class SraRetryConfigTest extends Specification {

    def 'should create retry config'() {

        expect:
        new SraRetryConfig().delay == Duration.of('500ms')
        new SraRetryConfig().maxDelay == Duration.of('30s')
        new SraRetryConfig().maxAttempts == 3
        new SraRetryConfig().jitter == 0.25d

        and:
        new SraRetryConfig([maxAttempts: 20]).maxAttempts == 20
        new SraRetryConfig([delay: '1s']).delay == Duration.of('1s')
        new SraRetryConfig([maxDelay: '1m']).maxDelay == Duration.of('1m')
        new SraRetryConfig([jitter: '0.5']).jitter == 0.5d

    }
}
