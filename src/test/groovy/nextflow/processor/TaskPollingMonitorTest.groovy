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

package nextflow.processor
import java.util.concurrent.ArrayBlockingQueue

import nextflow.Session
import nextflow.util.Duration
import nextflow.util.RateUnit
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskPollingMonitorTest extends Specification {

    def testCreate() {

        setup:
        def name = 'hello'
        def session = new Session( executor: [pollInterval: '1h', queueSize: 11, dumpInterval: '3h'] )

        def defSize = 99
        def defPollDuration = Duration.of('44s')
        when:
        def monitor = TaskPollingMonitor.create(session, name, defSize, defPollDuration)
        then:
        monitor.name == 'hello'
        monitor.dispatcher == session.dispatcher
        monitor.pollIntervalMillis == Duration.of('1h').toMillis()
        monitor.capacity == 11
        monitor.dumpInterval ==  Duration.of('3h')
        monitor.runningQueue instanceof ArrayBlockingQueue

    }

    def 'should create a rate limiter for the given rate format'() {

        given:
        def session = Mock(Session)
        def monitor = new TaskPollingMonitor(name:'local', session: session, pollInterval: '1s', capacity: 100)
        
        when:
        def limit = monitor.createSubmitRateLimit()
        then:
        1 * session.getExecConfigProp('local','submitRateLimit', null) >> RATE
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


    def 'check equals and hash code' () {
        expect:
        new RateUnit(2.1) == new RateUnit(2.1)
        new RateUnit(2.1).hashCode() == new RateUnit(2.1).hashCode()
        new RateUnit(2.1) != new RateUnit(3.3)
        new RateUnit(2.1).hashCode() != new RateUnit(3.3).hashCode()
    }

    def 'should stringify' () {
        expect:
        new RateUnit(0.1).toString() == '0.10/sec'
        new RateUnit(2.1).toString() == '2.10/sec'
        new RateUnit(123.4).toString() == '123.40/sec'
    }

}
