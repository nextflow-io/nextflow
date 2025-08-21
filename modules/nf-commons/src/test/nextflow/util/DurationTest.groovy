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

import java.time.Instant
import java.time.LocalDateTime
import java.time.temporal.ChronoUnit
import java.time.temporal.UnsupportedTemporalTypeException

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DurationTest extends Specification {

    long SEC = 1000
    long MIN = 60 * SEC
    long HOUR = 60 * MIN
    long DAY = 24 * HOUR

    def 'test create by string' () {

        expect:
        Duration.of('123 millis').toMillis() == 123
        Duration.of('123 ms').toMillis() == 123
        Duration.of('123ms').toMillis() == 123

        Duration.of('5 seconds').toSeconds() == 5
        Duration.of('5 second').toSeconds() == 5
        Duration.of('5 sec').toSeconds() == 5
        Duration.of('5 s').toSeconds() == 5
        Duration.of('5s').toSeconds() == 5

        Duration.of('5 minutes').toSeconds() == 300
        Duration.of('5 minute').toSeconds() == 300
        Duration.of('5 min').toSeconds() == 300
        Duration.of('5min').toSeconds() == 300

        Duration.of('5 hours').toMinutes() == 300
        Duration.of('5 hour').toMinutes() == 300
        Duration.of('5 h').toMinutes() == 300
        Duration.of('5h').toMinutes() == 300

        Duration.of('1 days').toHours() == 24
        Duration.of('1 day').toHours() == 24
        Duration.of('1 d').toHours() == 24
        Duration.of('1d').toHours() == 24
        Duration.of('1days').toHours() == 24
        Duration.of('1day').toHours() == 24
        Duration.of('1d').toHours() == 24

        Duration.of('30d').toMillis() == 30 * DAY

    }

    def 'should not be parsed as a duration' () {
        when:
        new Duration('live_in_3d')
        then:
        thrown(IllegalArgumentException)

        when:
        new Duration('/path/to/samples/2016-06-05_21:04:05/sample.bam')
        then:
        thrown(IllegalArgumentException)
    }

    def 'should parse multi unit time format'() {
        expect:
        Duration.of('1d 2h').toMillis() ==  DAY + 2 * HOUR
        Duration.of('1 d 2 h').toMillis() ==  DAY + 2 * HOUR
        Duration.of('2d3h4m').toMillis() ==  2 * DAY + 3 * HOUR + 4 * MIN
        Duration.of('2d 3h 4m').toMillis() ==  2 * DAY + 3 * HOUR + 4 * MIN
    }

    def 'should parse float time' () {
        expect:
        Duration.of('10.5 s').toMillis() == 10_500
        Duration.of('10.5 m').toSeconds() == 630
    }

    def 'should parse legacy time format string' () {

        expect:
        Duration.of('1:0:0').toString() == '1h'
        Duration.of('01:00:00').toString() == '1h'
        Duration.of('10:00:00').toString() == '10h'
        Duration.of('01:02:03').toString() == '1h 2m 3s'
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
        duration.toString() == '5m'

    }

    def 'test toString' () {

        expect:
        new Duration(100).toString() == '100ms'
        new Duration(1_000).toString() == '1s'
        new Duration(1_100).toString() == '1.1s'
        new Duration(32_300).toString() == '32.3s'
        new Duration(61 * 1000).toString() == '1m 1s'
        new Duration(61 * 1000 + 200).toString() == '1m 1s'
        new Duration(61 * 1000 + 800).toString() == '1m 2s'
        new Duration(60 * 60 * 1000 + 1000).toString() == '1h 1s'
        new Duration(25 * 60 * 60 * 1000 + 1000).toString() == '1d 1h 1s'

    }

    def 'test equals and compare' () {
        expect:
        new Duration('1h') == new Duration('1h')
        new Duration('10min') < new Duration('11min')
        new Duration('20min') > new Duration('11min')

    }

    def 'should not convert UUID number' () {

        when:
        new Duration('10833d95-1546-4BDF-AEA6-9F9676571854')
        then:
        thrown(IllegalArgumentException)

    }

    def 'should add time' () {

        expect:
        new Duration('1 hour') + new Duration('3 hours')  == new Duration('4 hours')

    }

    def 'should subtract time' () {

        expect:
        new Duration('4 hour') - new Duration('3 hours')  == new Duration('1 hour')

    }

    def 'should multiple time' () {
        expect:
        new Duration('2 hour') * 3 == new Duration('6 hours')
        new Duration('6 hours') * 1.5 == new Duration('9 hours')
        // `multiply` a number by a MemoryUnit is implemented by `NumberDelegatingMetaClass`
    }

    def 'should divide time' () {
        expect:
        new Duration('6 hour') / 2 == new Duration('3 hours')
        new Duration('9 hours') / 1.5 == new Duration('6 hours')
    }

    def 'should validate groovy truth' () {
        expect:
        !new Duration(0)
        new Duration(1)
    }


    def 'should validate duration between' () {

        given:
        def start = Instant.now()
        def end = start.plusMillis(1000)

        expect:
        Duration.between(start, end) == Duration.of('1sec')
    }

    // TemporalAmount interface tests

    def 'should get value for supported temporal units' () {
        given:
        def duration = new Duration('2d 3h 4m 5s 123ms')

        expect:
        duration.get(ChronoUnit.MILLIS) == duration.toMillis()
        duration.get(ChronoUnit.SECONDS) == duration.toSeconds()
        duration.get(ChronoUnit.MINUTES) == duration.toMinutes()
        duration.get(ChronoUnit.HOURS) == duration.toHours()
        duration.get(ChronoUnit.DAYS) == duration.toDays()
    }

    def 'should throw exception for unsupported temporal unit' () {
        given:
        def duration = new Duration('1h')

        when:
        duration.get(ChronoUnit.WEEKS)

        then:
        thrown(UnsupportedTemporalTypeException)
    }

    def 'should return list of supported units' () {
        given:
        def duration = new Duration('1h')

        when:
        def units = duration.getUnits()

        then:
        units.size() == 5
        units.contains(ChronoUnit.DAYS)
        units.contains(ChronoUnit.HOURS)
        units.contains(ChronoUnit.MINUTES)
        units.contains(ChronoUnit.SECONDS)
        units.contains(ChronoUnit.MILLIS)
        // Should be ordered from longest to shortest
        units[0] == ChronoUnit.DAYS
        units[4] == ChronoUnit.MILLIS
    }

    def 'should add duration to temporal' () {
        given:
        def duration = new Duration('2h 30m')
        def dateTime = LocalDateTime.of(2023, 1, 1, 10, 0, 0)

        when:
        def result = duration.addTo(dateTime)

        then:
        result == LocalDateTime.of(2023, 1, 1, 12, 30, 0)
    }

    def 'should subtract duration from temporal' () {
        given:
        def duration = new Duration('1h 15m')
        def dateTime = LocalDateTime.of(2023, 1, 1, 12, 30, 0)

        when:
        def result = duration.subtractFrom(dateTime)

        then:
        result == LocalDateTime.of(2023, 1, 1, 11, 15, 0)
    }

    def 'should add duration to instant' () {
        given:
        def duration = new Duration('5s')
        def instant = Instant.ofEpochMilli(1000)

        when:
        def result = duration.addTo(instant)

        then:
        result == Instant.ofEpochMilli(6000)
    }

    def 'should subtract duration from instant' () {
        given:
        def duration = new Duration('2s')
        def instant = Instant.ofEpochMilli(5000)

        when:
        def result = duration.subtractFrom(instant)

        then:
        result == Instant.ofEpochMilli(3000)
    }


}
