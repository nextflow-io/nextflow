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

package nextflow.trace

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WorkflowStatsTest extends Specification {


    def 'should return compute time string' () {

        given:
        WorkflowStats stats

        when:
        stats = new WorkflowStats(succeedMillis: 3_600_000)
        then:
        stats.getComputeTimeString() == '1.0'

        when:
        stats = new WorkflowStats(succeedMillis: 10_000_000)
        then:
        stats.getComputeTimeString() == '2.8'

        when:
        stats = new WorkflowStats(succeedMillis: 80_000_000, cachedMillis: 20_000_000)
        then:
        stats.getComputeTimeString() == '27.8 (20% cached)'

        when:
        stats = new WorkflowStats(succeedMillis: 80_000_000, failedMillis: 20_000_000)
        then:
        stats.getComputeTimeString() == '27.8 (20% failed)'

        when:
        stats = new WorkflowStats(succeedMillis: 80_000_000, failedMillis: 5_000_000, cachedMillis: 15_000_000)
        then:
        stats.getComputeTimeString() == '27.8 (15% cached, 5% failed)'

        when:
        stats = new WorkflowStats(succeedMillis: 180_000)
        then:
        stats.getComputeTimeString() == '0.1'

        when:
        stats = new WorkflowStats(succeedMillis: 120_000)
        then:
        stats.getComputeTimeString() == '(a few seconds)'

        when:
        stats = new WorkflowStats(succeedMillis: 120_000_000_000)
        then:
        stats.getComputeTimeString() == "33'333.3"
    }

    def 'should return cpu time' () {
        given:
        long seconds
        def stats = new WorkflowStats()

        when:
        seconds = stats.getCpuTime( new TraceRecord() )
        then:
        seconds == 0

        when:
        seconds = stats.getCpuTime( new TraceRecord([realtime: 35_000, cpus: 1]) )
        then:
        seconds == 35_000

        when:
        seconds = stats.getCpuTime( new TraceRecord([realtime: 10_000, cpus: 4]) )
        then:
        seconds == 40_000
    }

    def 'should update cached stats' () {

        given:
        def stats = new WorkflowStats()

        when:
        stats.updateTasksCached( new TraceRecord([realtime: 30_000, cpus: 1]) )
        stats.updateTasksCached( new TraceRecord([realtime: 10_000, cpus: 4]) )
        stats.updateTasksCached( new TraceRecord([realtime: 5_000, cpus: 1]) )
        then:
        stats.cachedCount == 3
        stats.cachedDuration.seconds == 75

    }

    def 'should update completed stats' () {
        given:
        def stats = new WorkflowStats()

        /*
         * should increment `succeed` stats
         */
        when:
        stats.updateTasksCompleted( new TraceRecord([status: 'COMPLETED', realtime: 30_000, cpus: 2]) )
        then:
        stats.succeedCount == 1
        stats.failedCount == 0
        stats.ignoredCount == 0
        stats.succeedDuration.seconds == 60
        stats.failedDuration.seconds == 0

        /*
         * should increment `succeed` stats
         */
        when:
        stats.updateTasksCompleted( new TraceRecord([status: 'COMPLETED', realtime: 10_000, cpus: 1]) )
        then:
        stats.succeedCount == 2
        stats.failedCount == 0
        stats.ignoredCount == 0
        stats.succeedDuration.seconds == 70
        stats.failedDuration.seconds == 0

        /*
         * should increment `failed` stats
         */
        when:
        stats.updateTasksCompleted( new TraceRecord([status: 'FAILED', realtime: 10_000, cpus: 1]) )
        then:
        stats.succeedCount == 2
        stats.failedCount == 1
        stats.ignoredCount == 0
        stats.succeedDuration.seconds == 70
        stats.failedDuration.seconds == 10

        /*
         * should increment `ignored` counter and failed time
         */
        when:
        stats.updateTasksCompleted( new TraceRecord([status: 'FAILED', realtime: 10_000, cpus: 1, error_action: 'IGNORE']) )
        then:
        stats.succeedCount == 2
        stats.failedCount == 1
        stats.ignoredCount == 1
        stats.succeedDuration.seconds == 70
        stats.failedDuration.seconds == 20
    }


    def 'should get cpu time' () {
        given:
        def stats = new WorkflowStats()
        def record = Mock(TraceRecord)
        long result

        when:
        result = stats.getCpuTime(record)
        then:
        1 * record.get('realtime') >> 1_000
        1 * record.get('cpus') >> 1
        result == 1000

        when:
        result = stats.getCpuTime(record)
        then:
        1 * record.get('realtime') >> 2_000
        1 * record.get('cpus') >> 4
        result == 8_000

        when:
        result = stats.getCpuTime(record)
        then:
        1 * record.get('realtime') >> 2_000
        1 * record.get('cpus') >> null
        result == 2_000

        when:
        result = stats.getCpuTime(record)
        then:
        1 * record.get('realtime') >> null
        1 * record.get('cpus') >> null
        result == 0
    }

}
