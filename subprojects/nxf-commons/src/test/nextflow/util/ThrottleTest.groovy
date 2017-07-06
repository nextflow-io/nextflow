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


}
