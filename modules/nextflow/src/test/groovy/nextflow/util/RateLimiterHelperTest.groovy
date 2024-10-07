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
 */
package nextflow.util

import spock.lang.Specification

/**
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class RateLimiterHelperTest extends Specification {

    def 'should create a rate limiter for the given rate format'() {

        when:
        def limit = RateLimiterHelper.createRateLimit(RATE)
        then:
        limit ? Math.round(limit.getRate()) : null == EXPECTED

        where:
        RATE            | EXPECTED
        '1'             | 1                         // 1 per second
        '5'             | 5                         // 5 per second
        '100 min'       | 100 / 60                  // 100 per minute
        '100 / 1 s'     | 100                       // 100 per second
        '100 / 2 s'     | 50                        // 100 per 2 seconds
        '200 / sec'     | 200                       // 200 per second
        '600 / 5'       | 600i / 5l as double       // 600 per 5 seconds
        '600 / 5min'    | 600 / (5 * 60)            // 600 per 5 minutes

    }
}
