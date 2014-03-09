package nextflow.script

import static test.TestParser.parse

import nextflow.executor.AbstractExecutor
import nextflow.executor.LocalExecutor
import nextflow.executor.SgeExecutor
import nextflow.processor.TaskProcessor
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
        BaseScript.loadExecutorClass('oge') == SgeExecutor

    }


    def testSupportType() {

        when:
        BaseScript.isTypeSupported(ScriptType.GROOVY, 'xxx')
        then:
        thrown(IllegalArgumentException)

        expect:
        // by default only script are supported
        BaseScript.isTypeSupported(ScriptType.SCRIPTLET, AbstractExecutor,)
        !BaseScript.isTypeSupported(ScriptType.GROOVY, AbstractExecutor)
        // same for grid
        BaseScript.isTypeSupported(ScriptType.SCRIPTLET, SgeExecutor)
        !BaseScript.isTypeSupported(ScriptType.GROOVY, SgeExecutor)

        // local supports both
        BaseScript.isTypeSupported(ScriptType.SCRIPTLET, LocalExecutor)
        BaseScript.isTypeSupported(ScriptType.GROOVY, LocalExecutor)

        // repeat for instances
        BaseScript.isTypeSupported(ScriptType.SCRIPTLET, new SgeExecutor() )
        !BaseScript.isTypeSupported(ScriptType.GROOVY, new SgeExecutor())
        BaseScript.isTypeSupported(ScriptType.SCRIPTLET, new LocalExecutor())
        BaseScript.isTypeSupported(ScriptType.GROOVY, new LocalExecutor())
    }

    def testVarRefs( ) {

        setup:
        def text = '''

        String x = 1

        @Field
        String = 'Ciao'

        z = 'str'

        process hola {

          /
          println $x + $y + $z
          /
        }
        '''
        when:
        TaskProcessor process = parse(text).run()
        then:
        process.taskBody.valRefs == [
                new TokenValRef('x', 13, 20),
                new TokenValRef('y', 13, 25),
                new TokenValRef('z', 13, 30) ] as Set

        process.taskBody.getValNames().sort() == ['x','y','z']
    }

}
