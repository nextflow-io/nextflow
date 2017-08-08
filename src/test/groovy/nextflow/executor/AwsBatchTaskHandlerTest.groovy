package nextflow.executor

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsBatchTaskHandlerTest extends Specification {

    def 'should normalise a task name' () {
        given:
        def handler = [:] as AwsBatchTaskHandler
        expect:
        handler.normalizeName('foo') == 'foo'
        handler.normalizeName('foo (12)') == 'foo_12'
    }

}
