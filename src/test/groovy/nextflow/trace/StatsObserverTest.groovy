package nextflow.trace
import nextflow.processor.TaskHandler
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class StatsObserverTest extends Specification {

    def 'should update tasks completed' () {
        given:
        def stats = Mock(WorkflowStats)
        def observer = new StatsObserver(stats: stats)

        def RECORD1 = Mock(TraceRecord)
        def RECORD2 = Mock(TraceRecord)

        def HANDLER1 = Mock(TaskHandler)
        def HANDLER2 = Mock(TaskHandler)

        when:
        observer.onProcessComplete(HANDLER1)
        then:
        1 * HANDLER1.getTraceRecord() >> RECORD1
        1 * stats.updateTasksCompleted(RECORD1) >> null

        when:
        observer.onProcessCached(HANDLER2)
        then:
        1 * HANDLER2.getTraceRecord() >> RECORD2
        1 * stats.updateTasksCached(RECORD2) >> null
    }


}
