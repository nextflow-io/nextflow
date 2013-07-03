package nextflow.script

import nextflow.executor.GenericGridExecutor
import nextflow.executor.LocalExecutor
import nextflow.executor.SgeExecutor
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BaseScriptTest extends Specification {


    def 'test loadStrategyClass' () {

        expect:
        BaseScript.loadStrategyClass(null) == LocalExecutor
        BaseScript.loadStrategyClass('local') == LocalExecutor
        BaseScript.loadStrategyClass('sge') == SgeExecutor
        BaseScript.loadStrategyClass( GenericGridExecutor.name ) == GenericGridExecutor

    }


}
