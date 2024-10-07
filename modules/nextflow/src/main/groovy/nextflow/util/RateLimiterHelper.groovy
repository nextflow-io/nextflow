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

import com.google.common.util.concurrent.RateLimiter
import groovy.transform.CompileStatic


@CompileStatic
class RateLimiterHelper {

    private static String RATE_FORMAT = ~/^(\d+\.?\d*)\s*([a-zA-Z]*)/

    static RateLimiter createRateLimit(String limit) {

        def tokens = limit.tokenize('/')
        if( tokens.size() == 2 ) {
            /*
             * the rate limit is provide num of task over a duration
             * - eg. 100 / 5 min
             * - ie. max 100 task per 5 minutes
             */

            final X = tokens[0].trim()
            final Y = tokens[1].trim()

            return newRateLimiter(X, Y, limit)
        }

        /*
         * the rate limit is provide as a duration
         * - eg. 200 min
         * - ie. max 200 task per minutes
         */

        final matcher = (limit =~ RATE_FORMAT)
        if( !matcher.matches() )
            throw new IllegalArgumentException("Invalid submit-rate-limit value: $limit -- It must be provide using the following format `num request sec|min|hour` eg. 10 sec ie. max 10 tasks per second")

        final num = matcher.group(1) ?: '_'
        final unit = matcher.group(2) ?: 'sec'

        return newRateLimiter(num, "1 $unit", limit)
    }

    private static RateLimiter newRateLimiter( String X, String Y, String limit ) {
        if( !X.isInteger() )
            throw new IllegalArgumentException("Invalid submit-rate-limit value: $limit -- It must be provided using the following format `num request / duration` eg. 10/1s")

        final num = Integer.parseInt(X)
        final duration = Y.isInteger() ? Duration.of( Y+'sec' ) : ( Y[0].isInteger() ? Duration.of(Y) : Duration.of('1'+Y) )
        long seconds = duration.toSeconds()
        if( !seconds )
            throw new IllegalArgumentException("Invalid submit-rate-limit value: $limit -- The interval must be at least 1 second")

        return RateLimiter.create( num / seconds as double )
    }
}
