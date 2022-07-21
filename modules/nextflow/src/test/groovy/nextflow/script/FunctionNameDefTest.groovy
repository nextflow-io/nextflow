package nextflow.script

import nextflow.NF
import spock.lang.Specification
/**
 *
 * @author Jorge Aguilera <jorge.aguilera@seqera.io>
 */
class FunctionNameDefTest extends Specification {

    static class TestScript extends BaseScript {

        def foo(String x='foo') {
            return "I'm $x"
        }

        def foo(Integer x) {
            return "I'm an integer $x"
        }

        def bar(String x) {
            return "I'm $x"
        }

        @Override
        protected Object runScript() { return null }
    }

    def setupSpec(){
        NF.init()
    }

    def 'should create a function def' () {
        given:
        def script = new TestScript()
        def f1 = TestScript.class.getMethod('foo')
        def f1_1 = TestScript.class.getMethod('foo', [String] as Class[])
        def f1_2 = TestScript.class.getMethod('foo', [Integer] as Class[])
        def f2 = TestScript.class.getMethod('bar', [String] as Class[])

        when:
        def func = new FunctionNameDef(script, 'foo')
        then:
        func.name == 'foo'
        func.invoke_o( new Object[]{} ) == "I'm foo"
        func.invoke_o( 'a string' ) == "I'm a string"
        func.invoke_o( 123 ) == "I'm an integer 123"

        when:
        func = new FunctionNameDef(script, 'bar')
        then:
        func.name == 'bar'
        func.invoke_o('Mr Bar') == "I'm Mr Bar"
    }


    def 'should clone with a new name' () {
        given:
        def script = new TestScript()
        def f1 = TestScript.class.getMethod('foo')

        when:
        def func = new FunctionDef(script, f1)
        def copy = func.cloneWithName('foo_alias')
        then:
        copy.name == 'foo_alias'
        copy.invoke_o() == "I'm foo"
    }

}
