/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.trace

import groovy.json.JsonSlurper
import spock.lang.Specification
import test.TestHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TraceRecordTest extends Specification {

    static long timestamp = 1412596982622

    def setupSpec() {
        TraceRecord.TIMEZONE = TimeZone.getTimeZone('UTC') // note: set the timezone to be sure the time string does not change on CI test servers
    }

    def 'test get field'() {

        given:
        def record = new TraceRecord()
        record.task_id = 'hola'
        record.native_id = null
        record.submit = timestamp
        record.duration = '2000'
        record.'%cpu' = '5.00'
        record.rss = '1024'
        record.queue = 'big'

        expect:
        record.getFmtStr(name, converter) == expected

        where:
        name        | converter | expected
        'task_id'   | null      | 'hola'
        'native_id' | null      | '-'
        'submit'    | null      | '2014-10-06 12:03:02.622'
        'submit'    | 'num'     | timestamp.toString()
        'duration'  | null      | '2s'
        '%cpu'      | null      | '5.0%'
        'rss'       | null      | '1 KB'
        'queue'     | null      | 'big'

    }

    def 'test fmt str' () {

        expect:
        TraceRecord.fmtString(str, fmt) == expect

        where:
        str     | fmt   | expect
        null    | null  | '-'
        '0'     | null  | '0'
        'hello' | null  | 'hello'
        '123'   | null  | '123'
        456     | null  | '456'
        0       | null  | '0'
        Integer.MAX_VALUE     | null  | '-'

    }


    def 'test fmt memory'() {
        expect:
        TraceRecord.fmtMemory(value, fmt) == expect

        where:
        value    | fmt  | expect
        null     | null | '-'
        2048     | null | '2 KB'
        1024000  | null | '1000 KB'
        '0'      | null | '0'
        '100'    | null | '100 B'
        '1024'   | null | '1 KB'
        ' 2048 ' | null | '2 KB'
        'abc'    | null | 'abc'
    }

    def 'test fmt date'() {
        expect:
        TraceRecord.fmtDate(str, fmt) == expect

        where:
        str                     | fmt   | expect
        null                    | null  | '-'
        1408714875000           | null  | '2014-08-22 13:41:15.000'
    }

    def 'test fmt time'() {
        expect:
        TraceRecord.fmtTime(str, fmt) == expect

        where:
        str                     | fmt   | expect
        0                       | null  | '0ms'
        100                     | null  | '100ms'
        2000                    | null  | '2s'
        3600 * 1000 * 3         | null  | '3h'
        3600 * 1000 * 3 + 5000  | null  | '3h 5s'

    }


    def 'test fmt number'() {
        expect:
        TraceRecord.fmtNumber(str, fmt) == expect

        where:
        str                     | fmt   | expect
        0                       | null  | '0'
        100                     | null  | '100'
        '333'                   | null  | '333'
        null                    | null  | '-'

    }


    def 'test fmt percent'() {
        expect:
        TraceRecord.fmtPercent(str, fmt) == expect

        where:
        str                     | fmt   | expect
        0                       | null  | '0.0%'
        1                       | null  | '1.0%'
        '0.199'                 | null  | '0.2%'
        '100.991'               | null  | '101.0%'
        'abc'                   | null  | '-'
    }


    static final long KB = 1024L

    def 'should parse a trace file and return a TraceRecord object'() {

        given:
        def file = TestHelper.createInMemTempFile('trace')
        file.text = '''\
            nextflow.trace/v2
            realtime=12021
            %cpu=997
            rchar=50838
            wchar=317
            syscr=120
            syscw=14
            read_bytes=0
            write_bytes=0
            %mem=9
            vmem=323104
            rss=146536
            peak_vmem=323252
            peak_rss=197136
            '''.stripIndent().leftTrim()

        when:
        def handler = [:] as TraceRecord
        def trace = handler.parseTraceFile(file)

        then:
        trace.realtime == 12021
        trace.'%cpu' == 99.7
        trace.rchar == 50838
        trace.wchar == 317
        trace.syscr == 120
        trace.syscw == 14
        trace.read_bytes == 0
        trace.write_bytes == 0
        trace.'%mem' == 0.9
        trace.vmem == 323104 * KB
        trace.rss == 146536 * KB
        trace.peak_vmem == 323252 * KB
        trace.peak_rss == 197136 * KB

        trace.getFmtStr('%mem') == '0.9%'
        trace.getFmtStr('vmem') == '315.5 MB'
        trace.getFmtStr('rss') == '143.1 MB'
        trace.getFmtStr('peak_vmem') == '315.7 MB'
        trace.getFmtStr('peak_rss') == '192.5 MB'
    }

    def 'should parse a legacy trace file and return a TraceRecord object'() {

        given:
        def file = TestHelper.createInMemTempFile('trace')
        file.text =  '''
        pid state %cpu %mem vmem rss peak_vmem peak_rss rchar wchar syscr syscw read_bytes write_bytes
        1 0 10 20 11084 1220 21084 2220 4790 12 11 1 20 30
        153
        '''
                .leftTrim()

        when:
        def handler = [:] as TraceRecord
        def trace = handler.parseTraceFile(file)
        then:
        trace.'%cpu' == 1.0
        trace.'%mem' == 2.0
        trace.rss == 1220 * KB
        trace.vmem == 11084 * KB
        trace.peak_rss == 2220 * KB
        trace.peak_vmem == 21084 * KB
        trace.rchar == 4790
        trace.wchar == 12
        trace.syscr == 11
        trace.syscw ==  1
        trace.read_bytes == 20
        trace.write_bytes == 30
        trace.realtime == 153

    }


    def 'should render record json' () {

        given:
        def record = new TraceRecord()
        record.task_id = 'hola'
        record.native_id = null
        record.submit = timestamp
        record.duration = '2000'
        record.'%cpu' = '5.00'
        record.rss = '1024'
        record.queue = 'big'
        record.cpus = 4
        record.time = 3_600_000L
        record.memory = 1024L * 1024L * 1024L * 8L

        when:
        def json = new JsonSlurper().parseText(record.renderJson().toString())

        then:
        json.task_id == 'hola'
        json.native_id == '-'
        json.submit == '2014-10-06 12:03:02.622'
        json.duration == '2s'
        json.'%cpu' == '5.0%'
        json.rss == '1 KB'
        json.queue == 'big'
        json.cpus == '4'
        json.time == '1h'
        json.memory == '8 GB'

    }


    def 'should remove secret key' () {
        given:
        def rec = new TraceRecord()
        expect:
        rec.secureEnvString('a=b') == 'a=b'
        and:
        rec.secureEnvString('aws_key=12345') == 'aws_key=[secure]'
        rec.secureEnvString('AWS_KEY=12345') == 'AWS_KEY=[secure]'

        rec.secureEnvString('''\
                foo=hello    
                aws_key=d7sds89
                git_token=909s-ds-'''
                .stripIndent() ) ==
                '''\
                foo=hello    
                aws_key=[secure]
                git_token=[secure]'''.stripIndent()

    }

    def 'should store safe env' () {
        given:
        def rec = new TraceRecord()
        when:
        rec.env = 'aws_key=1234'
        then:
        rec.store.env == 'aws_key=[secure]'
    }

    def 'should retrieve safe env' () {
        given:
        def rec = new TraceRecord()
        when:
        rec.store.env = 'aws_key=1234'
        then:
        rec.env == 'aws_key=[secure]'
    }

}
