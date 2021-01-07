package nextflow.cloud.azure.batch

import nextflow.cloud.azure.config.AzConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzBatchServiceTest extends Specification {

    def 'should make job id'() {
        given:
        def task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) {
                getName() >> NAME
            }
        }
        and:
        def exec = Mock(AzBatchExecutor) {
            getConfig() >> new AzConfig([:])
        }
        and:
        def svc = new AzBatchService(exec)

        expect:
        svc.makeJobId(task) == EXPECTED

        where:
        NAME        | EXPECTED
        'foo'       | 'nf-job-foo'
        'foo  bar'  | 'nf-job-foo_bar'
    }
}
