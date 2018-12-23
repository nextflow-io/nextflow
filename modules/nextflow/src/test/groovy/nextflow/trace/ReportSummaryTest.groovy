/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ReportSummaryTest extends Specification {

    private double round ( double value ) {
        Math.round(value * 100 ) / 100
    }

    private quantile( List items, int q0 ) {
        def stats = new ReportSummary.Summary({ it as Double})
        round(stats.quantile(items.findAll{ it!=null } .sort(), q0))
    }

    private mean( List items ) {
        def notNull = items.findAll{ it!=null }
        round(notNull.sum() / notNull.size())
    }

    def 'should compute summary data' () {
        given:
        def summary = new ReportSummary.Summary({ it.get('%cpu') as Double })

        def values = [10,77,33,98,12,87,23,12,76,21]
        values.eachWithIndex { v, i ->
            def t = new TraceRecord('%cpu': v, name: "task-$i")
            summary.add(t)
        }

        when:
        def result = summary.compute()
        then:
        result == [
                mean: 44.9,
                min:10,
                q1: 14.25,
                q2: 28.0,
                q3: 76.75,
                max: 98,
                minLabel: 'task-0',
                maxLabel: 'task-3',
                q1Label: 'task-7',
                q2Label: 'task-6',
                q3Label: 'task-8'
        ]
    }

    def 'should compute empty summary' () {
        given:
        def summary = new ReportSummary.Summary({ it as Double })
        when:
        def result = summary.compute()
        then:
        result == null
    }

    def 'should get report summary stats' () {
        given:
        def CPU = [10,77,33,98,12,87,23,12,76,21]
        def MEM = [65,34,23,65,11,54,78,87,32,21]
        def TIME = [22,21,76,78,07,32,99,32,11,12]
        def READS = [87,32,65,87,23,11,19,87,29,33]
        def WRITES = [83,88,99,32,65,null,32,null,87,43]
        def data = [
            '%cpu': CPU,
            'vmem': MEM,
            realtime: TIME,
            'read_bytes': READS,
            'write_bytes': WRITES,
        ]

        def report = new ReportSummary()

        10.times { int index ->
            def t=new TraceRecord([process:'foo', name:"foo-$index"]);
            data.each { name, col ->  t.put(name, col[index])}
            // add a new record
            report.add(t)
        }


        when:
        def cpu = report.compute('cpu')
        def mem = report.compute('mem')
        def time = report.compute('time')
        def reads = report.compute('reads')
        def writes = report.compute('writes')
        then:

        cpu.min == quantile(CPU, 0)
        cpu.q1 == quantile(CPU, 25)
        cpu.q2 == quantile(CPU, 50)
        cpu.q3 == quantile(CPU, 75)
        cpu.max == quantile(CPU, 100)
        cpu.mean == mean(CPU)
        cpu.minLabel == 'foo-0'
        cpu.maxLabel == 'foo-3'
        cpu.q1Label == 'foo-7'
        cpu.q2Label == 'foo-6'
        cpu.q3Label == 'foo-8'

        mem.min == quantile(MEM, 0)
        mem.q1 == quantile(MEM, 25)
        mem.q2 == quantile(MEM, 50)
        mem.q3 == quantile(MEM, 75)
        mem.max == quantile(MEM, 100)
        mem.mean == mean(MEM)

        time.min == quantile(TIME, 0)
        time.q1 == quantile(TIME, 25)
        time.q2 == quantile(TIME, 50)
        time.q3 == quantile(TIME, 75)
        time.max == quantile(TIME, 100)
        time.mean == mean(TIME)

        reads.min == quantile(READS, 0)
        reads.q1 == quantile(READS, 25)
        reads.q2 == quantile(READS, 50)
        reads.q3 == quantile(READS, 75)
        reads.max == quantile(READS, 100)
        reads.mean == mean(READS)

        writes.min == quantile(WRITES, 0)
        writes.q1 == quantile(WRITES, 25)
        writes.q2 == quantile(WRITES, 50)
        writes.q3 == quantile(WRITES, 75)
        writes.max == quantile(WRITES, 100)
        writes.mean == mean(WRITES)
    }


    @Unroll
    def 'should calc quantiles' () {

        given:
        def stats = new ReportSummary.Summary({ it as Double })

        expect:
        stats.quantiles(NUMBERS.sort()) == EXPECTED

        where:
        EXPECTED                        | NUMBERS
        [1, 3.5, 6, 8.5, 11]            |  [1 , 3 , 5 , 7 , 9, 11]
        [7, 20.25, 37.50, 39.75, 41]    |  [7, 15, 36, 39, 40, 41 ]
        [10, 12, 15, 19.50, 26]         |  [10, 11, 13, 15, 16, 23, 26]
        [3, 7, 12, 14, 21 ]             |  [3, 5, 7, 8, 12, 13, 14, 18, 21]
        [6, 25.50, 40, 42.50, 49]       |  [6, 7, 15, 36, 39, 40, 41, 42, 43, 47, 49]
        [7, 20.25, 37.50, 39.75, 41]    |  [7, 15, 36, 39, 40, 41]
        [10, 12.50, 15.50, 23.75, 32]   |  [10, 11, 13, 15, 16, 23, 26, 32]
        [11, 244, 605, 787, 998]        |  [11,99,104,107,112,153,176,177,180,185,188,189,244,255,296,312,340,361,479,531,536,560,568,581,605,642,654,671,693,700,709,710,728,777,778,783,787,790,795,806,812,816,832,861,873,878,905,939,998]
        [7, 202.75, 593, 786, 998]      |  [7,11,99,104,107,112,153,176,177,180,185,188,189,244,255,296,312,340,361,479,531,536,560,568,581,605,642,654,671,693,700,709,710,728,777,778,783,787,790,795,806,812,816,832,861,873,878,905,939,998]

    }

    def 'should set q1,q2,q3 labels' () {
        given:
        def records = []
        [3, 5, 7, 8, 12, 13, 14, 18, 21].eachWithIndex { v, i ->
            records << new TraceRecord([realtime: v, name: "task-${i+1}"])
        }

        def mapping = { TraceRecord r -> r.get('realtime') as Double }
        def stats = new ReportSummary.Summary(mapping)
        records.sort(mapping)

        expect:
        stats.quantile(records, 0) == 3
        stats.quantile(records, 25) == 7
        stats.quantile(records, 50) == 12
        stats.quantile(records, 75) == 14
        stats.quantile(records, 100) == 21

        stats.minLabel == 'task-1'
        stats.q1Label == 'task-3'
        stats.q2Label == 'task-5'
        stats.q3Label == 'task-7'
        stats.maxLabel == 'task-9'

    }



}
