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
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ThrottlingExecutorTest extends Specification {

    def 'should increase rate limit' () {

        given:
        def change=0
        def opts = new ThrottlingExecutor.Options() .onRateLimitChange { change++ } .normalise()
        ThrottlingExecutor exec = Spy(ThrottlingExecutor, constructorArgs: [opts])
        def limiter = RateLimiter.create(1.0)

        when:
        exec.speedUp0(null)
        then:
        noExceptionThrown()
        limiter.rate == 1.0d

        when:
        exec.speedUp0(limiter)
        then:
        limiter.rate as float == 1.2d as float
        change == 1

    }

    def 'should NOT increase rate limit' () {

        given:
        def change=0
        def opts = new ThrottlingExecutor.Options()
                    .onRateLimitChange { change++ }
                    .withRampUp(5, 2)
                    .normalise()

        ThrottlingExecutor exec = Spy(ThrottlingExecutor, constructorArgs: [opts])
        def limiter = RateLimiter.create(4)

        when:
        exec.speedUp0(limiter)
        then:
        limiter.rate == 5 // <- use the max rate define in the options
        change == 1
    }


    def 'should decrease submit rate' () {
        given:
        def change=0
        def opts = new ThrottlingExecutor.Options()
                        .onRateLimitChange { change++ }
                        .normalise()
        ThrottlingExecutor exec = Spy(ThrottlingExecutor, constructorArgs: [opts])
        def limiter = RateLimiter.create(1.0)

        when:
        exec.backOff0(null)
        then:
        limiter.rate == 1.0d
        change == 0

        when:
        exec.backOff0(limiter)
        then:
        limiter.rate == 0.5
        change == 1 
    }

    def 'should not decrease submit rate' () {
        given:
        def change=0
        def opts = new ThrottlingExecutor.Options()
                        .onRateLimitChange { change++ }
                        .withBackOff(0.75)
                        .normalise()
        ThrottlingExecutor exec = Spy(ThrottlingExecutor, constructorArgs: [opts])
        def limiter = RateLimiter.create(1.0)

        when:
        exec.backOff0(limiter)
        then:
        limiter.rate == 0.75    // <-- use the min rate define in the setting
        change == 1

    }

}
