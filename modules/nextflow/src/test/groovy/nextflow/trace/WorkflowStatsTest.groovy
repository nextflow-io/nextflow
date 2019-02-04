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
        stats.getComputeTimeFmt() == '1.0'

        when:
        stats = new WorkflowStats(succeedMillis: 10_000_000)
        then:
        stats.getComputeTimeFmt() == '2.8'

        when:
        stats = new WorkflowStats(succeedMillis: 80_000_000, cachedMillis: 20_000_000)
        then:
        stats.getComputeTimeFmt() == '27.8 (20% cached)'

        when:
        stats = new WorkflowStats(succeedMillis: 80_000_000, failedMillis: 20_000_000)
        then:
        stats.getComputeTimeFmt() == '27.8 (20% failed)'

        when:
        stats = new WorkflowStats(succeedMillis: 80_000_000, failedMillis: 5_000_000, cachedMillis: 15_000_000)
        then:
        stats.getComputeTimeFmt() == '27.8 (15% cached, 5% failed)'

        when:
        stats = new WorkflowStats(succeedMillis: 180_000)
        then:
        stats.getComputeTimeFmt() == '0.1'

        when:
        stats = new WorkflowStats(succeedMillis: 120_000)
        then:
        stats.getComputeTimeFmt() == '(a few seconds)'

        when:
        stats = new WorkflowStats(succeedMillis: 120_000_000_000)
        then:
        stats.getComputeTimeFmt() == "33'333.3"
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

    def 'should return formatted counts' () {

        given:
        def stats = new WorkflowStats()

        when:
        stats.succeedCount = 12345
        then:
        stats.succeedCount == 12345
        stats.succeedCountFmt == "12'345"

        when:
        stats.cachedCount = 88774411
        then:
        stats.cachedCount == 88774411
        stats.cachedCountFmt == "88'774'411"

        when:
        stats.ignoredCount = 66332211
        then:
        stats.ignoredCount == 66332211
        stats.ignoredCountFmt == "66'332'211"

        when:
        stats.failedCount = 33776644
        then:
        stats.failedCount == 33776644
        stats.failedCountFmt == "33'776'644"
    }


    def 'should return task percents' () {
        given:
        def stats = new WorkflowStats(succeedCount: 20, cachedCount: 40, ignoredCount: 60, failedCount: 80)

        expect: 
        stats.getSucceedPct() == 10.0f
        stats.getCachedPct() == 20.0f
        stats.getIgnoredPct() == 30.0f
        stats.getFailedPct() == 40.0f

    }

}
