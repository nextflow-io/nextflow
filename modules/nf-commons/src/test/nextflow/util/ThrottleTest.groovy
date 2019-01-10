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
