package nextflow.script

import nextflow.executor.AbstractGridExecutor
import nextflow.executor.LocalExecutor
import nextflow.executor.SgeExecutor
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BaseScriptTest extends Specification {


    def 'test loadExecutor' () {

        expect:
        BaseScript.loadExecutorClass(null) == LocalExecutor
        BaseScript.loadExecutorClass('local') == LocalExecutor
        BaseScript.loadExecutorClass('sge') == SgeExecutor
        BaseScript.loadExecutorClass( AbstractGridExecutor.name ) == AbstractGridExecutor

    }


}
