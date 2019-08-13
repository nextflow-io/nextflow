package nextflow.script

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FunctionDefTest extends Specification {

    static class TestScript extends BaseScript {

        def foo() {
            return "I'm foo"
        }

        def bar(String x) {
            return "I'm $x"
        }

        @Override
        protected Object runScript() { return null }
    }


    def 'should create a function def' () {
        given:
        def script = new TestScript()
        def f1 = TestScript.class.getMethod('foo')
        def f2 = TestScript.class.getMethod('bar', [String] as Class[])

        when:
        def func = new FunctionDef(script, f1)
        then:
        func.name == 'foo'
        func.invoke_o() == "I'm foo"

        when:
        func = new FunctionDef(script, f2)
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
