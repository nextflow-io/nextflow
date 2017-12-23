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

}
