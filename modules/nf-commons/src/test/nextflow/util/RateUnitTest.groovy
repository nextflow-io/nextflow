/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class RateUnitTest extends Specification {


    @Unroll
    def 'should create rate value for string=#RATE' () {

        when:
        def limit = new RateUnit(RATE)
        then:
        limit.getRate() == EXPECTED

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

    def 'should create rate limit with a string' () {
        expect:
        RateUnit.of('1 /sec') == new RateUnit(1.0)
        RateUnit.of('1 /min') == new RateUnit( 1 / 60 )
    }


    def 'should create with a double value' () {
        expect:
        new RateUnit(5.32d).rate == 5.32d
    }

    def 'should create rate limiter' () {
        expect:
        new RateUnit(5.32d).rateLimiter.rate == 5.32d
    }
}
