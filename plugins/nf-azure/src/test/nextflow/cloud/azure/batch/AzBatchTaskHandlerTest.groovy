package nextflow.cloud.azure.batch

import nextflow.cloud.types.CloudMachineInfo
import nextflow.cloud.types.PriceModel
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.BashWrapperBuilder
import nextflow.executor.Executor
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.script.BaseScript
import nextflow.script.ProcessConfig
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


    def 'should create the trace record' () {
        given:
        def exec = Mock(Executor) { getName() >> 'azurebatch' }
        def processor = Mock(TaskProcessor)
        processor.getExecutor() >> exec
        processor.getName() >> 'foo'
        processor.getConfig() >> new ProcessConfig(Mock(BaseScript))
        def task = Mock(TaskRun)
        task.getProcessor() >> processor
        task.getConfig() >> GroovyMock(TaskConfig)
        def handler = Spy(AzBatchTaskHandler)
        handler.task = task
        handler.@taskKey = new AzTaskKey('job-123', 'nf-456')

        when:
        def trace = handler.getTraceRecord()
        then:
        1 * handler.isCompleted() >> false
        1 * handler.getMachineInfo() >> new CloudMachineInfo('Standard1', 'west-eu', PriceModel.standard)

        and:
        trace.native_id == 'job-123/nf-456'
        trace.executorName == 'azurebatch'
        trace.machineInfo.type == 'Standard1'
        trace.machineInfo.zone == 'west-eu'
        trace.machineInfo.priceModel == PriceModel.standard
    }

}
