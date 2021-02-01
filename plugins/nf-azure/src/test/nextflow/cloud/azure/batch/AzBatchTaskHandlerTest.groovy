package nextflow.cloud.azure.batch

import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.BashWrapperBuilder
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzBatchTaskHandlerTest extends Specification {

    def 'should validate config' () {
        when:
        def task = Mock(TaskRun) { getName() >> 'foo'; }
        and:
        new AzBatchTaskHandler(task: task)
                .validateConfiguration()
        then:
        def e = thrown(ProcessUnrecoverableException)
        e.message.startsWith('No container image specified for process foo')


        when:
        task = Mock(TaskRun) { getName() >> 'foo'; getContainer() >> 'ubuntu' }
        and:
        new AzBatchTaskHandler(task: task)
                .validateConfiguration()
        then:
        noExceptionThrown()
    }

    def 'should submit task' () {
        given:
        def builder = Mock(BashWrapperBuilder)
        def task = Mock(TaskRun)
        def azure = Mock(AzBatchService)
        and:
        def handler = Spy(AzBatchTaskHandler) {getBatchService() >> azure }
        handler.task = task
        when:
        handler.submit()
        then:
        1 * handler.createBashWrapper() >> builder
        1 * builder.build() >> null
        1 * azure.submitTask(task) >> null
        and:
        handler.getStatus() == TaskStatus.SUBMITTED
    }

}
