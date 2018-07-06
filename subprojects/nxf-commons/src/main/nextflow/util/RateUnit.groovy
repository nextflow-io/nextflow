/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.util

import com.google.common.util.concurrent.RateLimiter
import groovy.transform.EqualsAndHashCode

/**
 * Model a rate limit measure unit
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@EqualsAndHashCode(includeFields=true,includes='rate')
class RateUnit {

    private static String RATE_FORMAT = ~/^(\d+\.?\d*)\s*([a-zA-Z]*)/

    private double rate

    RateUnit(double rate) {
        this.rate = rate
    }

    RateUnit(String str) {
        this.rate = parse(str)
    }

    double getRate() { rate }

    RateLimiter getRateLimiter() { RateLimiter.create(rate) }

    protected double parse(String limit) {

        def tokens = limit.tokenize('/')
        if( tokens.size() == 2 ) {
            /*
             * the rate limit is provide num of task over a duration
             * - eg. 100 / 5 min
             * - ie. max 100 task per 5 minutes
             */

            final X = tokens[0].trim()
            final Y = tokens[1].trim()

            return parse0(X, Y, limit)
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

        return parse0(num, "1 $unit", limit)
    }


    private double parse0(String X, String Y, String limit ) {
        if( !X.isInteger() )
            throw new IllegalArgumentException("Invalid submit-rate-limit value: $limit -- It must be provide using the following format `num request / duration` eg. 10/1s")

        final num = Integer.parseInt(X)
        final duration = Y.isInteger() ? Duration.of( Y+'sec' ) : ( Y[0].isInteger() ? Duration.of(Y) : Duration.of('1'+Y) )
        long seconds = duration.toSeconds()
        if( !seconds )
            throw new IllegalArgumentException("Invalid submit-rate-limit value: $limit -- The interval must be at least 1 second")

        num / seconds as double
    }

    String toString() {
        String.format("%.2f", rate) + '/sec'
    }
}
