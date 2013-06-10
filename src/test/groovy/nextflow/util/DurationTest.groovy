/*
 * Copyright (c) 2012, the authors.
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
class DurationTest extends Specification {

    def 'test create by string' () {

        expect:
        Duration.create('123 millis').toMillis() == 123
        Duration.create('123 ms').toMillis() == 123
        Duration.create('123ms').toMillis() == 123

        Duration.create('5 seconds').toSeconds() == 5
        Duration.create('5 second').toSeconds() == 5
        Duration.create('5 sec').toSeconds() == 5
        Duration.create('5 s').toSeconds() == 5
        Duration.create('5s').toSeconds() == 5

        Duration.create('5 minutes').toSeconds() == 300
        Duration.create('5 minute').toSeconds() == 300
        Duration.create('5 min').toSeconds() == 300
        Duration.create('5min').toSeconds() == 300

        Duration.create('5 hours').toMinutes() == 300
        Duration.create('5 hour').toMinutes() == 300
        Duration.create('5 h').toMinutes() == 300
        Duration.create('5h').toMinutes() == 300

        Duration.create('1 days').toHours() == 24
        Duration.create('1 day').toHours() == 24
        Duration.create('1 d').toHours() == 24
        Duration.create('1d').toHours() == 24
        Duration.create('1days').toHours() == 24
        Duration.create('1day').toHours() == 24
        Duration.create('1d').toHours() == 24

    }

    def 'test format' () {

        when:
        def duration = new Duration('5min')

        then:
        duration.durationInMillis == 5 * 60 * 1000
        duration.toMillis() == 5 * 60 * 1000
        duration.toSeconds() == 5 * 60
        duration.toMinutes() == 5

        duration.format('ss') == '300'
        duration.format('mm:ss') == '05:00'
        duration.toString() == '5min'

    }

    def 'test toString' () {

        expect:
        new Duration(1000).toString() == '1sec'
        new Duration(61 * 1000).toString() == '1min 1sec'
        new Duration(60 * 60 * 1000 + 1000).toString() == '1hour 1sec'
        new Duration(25 * 60 * 60 * 1000 + 1000).toString() == '1day 1hour 1sec'

    }

    def 'test equals and compare' () {
        expect:
        new Duration('1h') == new Duration('1h')
        new Duration('10min') < new Duration('11min')
        new Duration('20min') > new Duration('11min')

    }

}
