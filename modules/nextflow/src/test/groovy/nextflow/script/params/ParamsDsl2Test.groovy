package nextflow.script.params

import nextflow.script.ScriptMeta
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.*
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class ParamsDsl2Test extends Dsl2Spec {

    def 'should allow unqualified stdin and stdout' () {

        given:
        def SCRIPT = '''

        process alpha {
            input:
            stdin
            output:
            stdout

            script:
            /echo foo/
        }

        workflow {}
        '''

        when:
        def script = loadScript(SCRIPT)
        and:
        def process = ScriptMeta.get(script).getProcess('alpha')

        then:
        def inputs = process.processConfig.getInputs()
        def outputs = process.processConfig.getOutputs()
        and:
        inputs.size() == 1
        inputs[0] instanceof StdInParam
        and:
        outputs.size() == 1
        outputs[0] instanceof StdOutParam

    }

    def 'should allow unqualified tuple stdin and stdout' () {

        given:
        def SCRIPT = '''

        process beta {
            input:
            tuple stdin, val(x)
            output:
            tuple stdout, path('z')

            script:
            /echo foo/
        }

        workflow {}
        '''

        when:
        def script = loadScript(SCRIPT)
        and:
        def process = ScriptMeta.get(script).getProcess('beta')

        then:
        def inputs = process.processConfig.getInputs()
        def outputs = process.processConfig.getOutputs()
        and:
        inputs.size() == 1
        and:
        def in_tuple = (TupleInParam) inputs.get(0)
        and:
        in_tuple.inner.size() == 2
        in_tuple.inner[0] instanceof StdInParam
        in_tuple.inner[1] instanceof ValueInParam

        and:
        def out_tuple = (TupleOutParam) outputs.get(0)
        and:
        out_tuple.inner.size() == 2
        out_tuple.inner[0] instanceof StdOutParam
        out_tuple.inner[1] instanceof FileOutParam

    }
}
