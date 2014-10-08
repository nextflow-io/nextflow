/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.trace

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TraceRecordTest extends Specification {

    static long timestamp = 1412596982622

    def testGet() {

        given:
        TraceRecord.getDateFormat().setTimeZone(TimeZone.getTimeZone('UTC')) // note: set the timezone to be sure the time string does not change on CI test servers§

        def record = new TraceRecord()
        record.task_id = 'hola'
        record.native_id = null
        record.submit = timestamp
        record.wall_time = '2000'
        record.'%cpu' = '5.00'
        record.rss = '1024'

        expect:
        record.get(name, converter) == expected

        where:
        name        | converter | expected
        'task_id'   | null      | 'hola'
        'native_id' | null      | '-'
        'submit'    | null      | '2014-10-06 12:03:02.622'
        'submit'    | 'num'     | timestamp.toString()
        'wall_time' | null      | '2s'
        '%cpu'     | null      | '5.00'
        'rss'       | null      | '1 KB'

    }

    @Unroll
    def testBytes() {
        expect:
        TraceRecord.fmtToMem(str, fmt) == expect

        where:
        str     | fmt   | expect
        null    | null  | '-'
        '0'     | null  | '0'
        '100'   | null  | '100 B'
        '1024'  | null  | '1 KB'
        ' 2048 '| null  | '2 KB'
        'abc'   | null  | 'abc'
    }

    def testDate() {
        expect:
        TraceRecord.fmtToDate(str, fmt) == expect

        where:
        str                     | fmt   | expect
        null                    | null  | '-'
        1408714875000           | null  | '2014-08-22 13:41:15.000'
    }

    def testTime() {
        expect:
        TraceRecord.fmtToTime(str, fmt) == expect

        where:
        str                     | fmt   | expect
        0                       | null  | '0ms'
        100                     | null  | '100ms'
        2000                    | null  |  '2s'
        3600 * 1000 * 3         | null  | '3h'
        3600 * 1000 * 3 + 5000  | null  | '3h 5s'

    }



//    def testPercent() {
//        expect:
//        TraceFileObserver.percent(0) == '0.00'
//        TraceFileObserver.percent(1) == '1.00'
//        TraceFileObserver.percent('0.123') == '0.12'
//        TraceFileObserver.percent('100.991') == '100.99'
//        TraceFileObserver.percent('abc') == '-'
//    }

}
