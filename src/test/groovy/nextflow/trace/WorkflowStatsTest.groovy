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
        stats = new WorkflowStats()
        stats.timeSucceed = 3600
        then:
        stats.getComputeTime() == '1.0'

        when:
        stats = new WorkflowStats()
        stats.timeSucceed = 10000
        then:
        stats.getComputeTime() == '2.8'

        when:
        stats = new WorkflowStats()
        stats.timeSucceed = 80_000
        stats.timeCached = 20_000
        then:
        stats.getComputeTime() == '27.8 (20% cached)'

        when:
        stats = new WorkflowStats()
        stats.timeSucceed = 80_000
        stats.timeFailed = 20_000
        then:
        stats.getComputeTime() == '27.8 (20% failed)'

        when:
        stats = new WorkflowStats()
        stats.timeSucceed = 80_000
        stats.timeFailed = 5_000
        stats.timeCached = 15_000
        then:
        stats.getComputeTime() == '27.8 (15% cached, 5% failed)'

    }

}
