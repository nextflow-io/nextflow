package nextflow.script

import nextflow.NF
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

    static class TestScriptOverloading extends BaseScript {

        def foo(String x='default') {
            return "I'm $x"
        }

        def foo(Integer i){
            return "I'm the integer $i"
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
        when:
        def func = new FunctionDef(script, 'foo')
        then:
        func.name == 'foo'
        func.invoke_o(new Object[]{}) == "I'm foo"

        when:
        func = new FunctionDef(script, 'bar')
        then:
        func.name == 'bar'
        func.invoke_o('Mr Bar') == "I'm Mr Bar"
    }


    def 'should clone with a new name' () {
        given:
        def script = new TestScript()

        when:
        def func = new FunctionDef(script, 'foo')
        def copy = func.cloneWithName('foo_alias')
        then:
        copy.name == 'foo_alias'
        copy.invoke_o(new Object[]{}) == "I'm foo"
    }

    def 'should create an overloaded function def' () {
        given:
        def script = new TestScriptOverloading()

        when:
        def func = new FunctionDef(script, "foo")
        then:
        func.name == 'foo'
        func.invoke_o(new Object[]{}) == "I'm default"

        and:
        func.invoke_o('foo') == "I'm foo"

        and:
        func.invoke_o(1) == "I'm the integer 1"

        when:
        func = new FunctionDef(script, 'bar')
        then:
        func.name == 'bar'
        func.invoke_o('Mr Bar') == "I'm Mr Bar"
    }
}
