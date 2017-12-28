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

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ThrottleTest extends Specification {

    def testThrottle() {

        when:
        int i = 0
        int result = -1
        5.times {
            result = Throttle.every('100ms') { (++i)*2 }
        }

        // it is invoked just one time
        // the result it is '2'
        then:
        i == 1
        result == 2

        when:
        i = 0
        result = -1
        5.times {
            result = Throttle.every('100ms') { (++i)*2 }
            sleep 30
        }

        // having a sleep of 20 millis and a delay of 100 millis, in 5 times it is called two times
        then:
        i == 2
        result == 4
    }

    def testThrottleWithSeed() {

        setup:

        when:
        int i = 0
        int result = -1
        5.times {
            result = Throttle.after('100ms', 33) { (++i)*2 }
        }
        then:
        // it is never invoked
        i == 0
        // and the 'seed' value it is returned
        result == 33


        when:
        i = 0
        result = -1
        5.times {
            result = Throttle.after('100ms') { (++i)*2 }
            sleep 30
        }

        // it is never invoked ust one time
        then:
        i == 1
        result == 2

    }

    def testCache() {
        given:
        int i=0
        int p=0
        int q=0

        when:
        p = Throttle.cache('X', 100) { ++i }
        p = Throttle.cache('X', 100) { ++i }
        p = Throttle.cache('X', 100) { ++i }
        then:
        p==1

        when:
        q = Throttle.cache('Y', 100) { ++i }
        q = Throttle.cache('Y', 100) { ++i }
        then:
        q == 2

        when:
        sleep 200
        q = Throttle.cache('Y', 100) { ++i }
        then:
        q == 3

    }

}
