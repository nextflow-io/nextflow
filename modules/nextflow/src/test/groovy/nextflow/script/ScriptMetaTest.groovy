package nextflow.script

import groovy.transform.InheritConstructors
import nextflow.exception.DuplicateModuleIncludeException
import test.Dsl2Spec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ScriptMetaTest extends Dsl2Spec {

    @InheritConstructors
    static class FooScript extends BaseScript {
        @Override
        protected Object runScript() { null }
    }

    def 'should return all defined names' () {

        given:
        def script = new FooScript(new ScriptBinding())

        def proc1 = new ProcessDef(script, Mock(Closure), 'proc1')
        def proc2 = new ProcessDef(script, Mock(Closure), 'proc2')
        def func1 = new FunctionDef(name: 'func1')
        def work1 = new WorkflowDef(name:'work1')

        def meta = new ScriptMeta(script)

        when:
        meta.addDefinition(func1)
        meta.addDefinition(work1)
        meta.addDefinition(proc1)
        meta.addDefinition(proc2)

        then:
        meta.getComponent('work1') == work1
        meta.getComponent('func1') == func1
        meta.getComponent('proc1') == proc1
        meta.getComponent('proc2') == proc2
        meta.getComponent('xxx') == null
        meta.getComponent('yyy') == null

        then:
        meta.getProcessNames() as Set == ['proc1','proc2'] as Set

    }

    def 'should add imports' () {

        given:
        def script1 = new FooScript(new ScriptBinding())
        def script2 = new FooScript(new ScriptBinding())
        def script3 = new FooScript(new ScriptBinding())
        def meta1 = new ScriptMeta(script1)
        def meta2 = new ScriptMeta(script2)
        def meta3 = new ScriptMeta(script3)

        // defs in the root script
        def func1 = new FunctionDef(name: 'func1')
        def proc1 = new ProcessDef(script1, Mock(Closure), 'proc1')
        def work1 = new WorkflowDef(name:'work1')
        meta1.addDefinition(proc1, func1, work1)

        // defs in the second script imported in the root namespace
        def func2 = new FunctionDef(name: 'func2')
        def proc2 = new ProcessDef(script2, Mock(Closure), 'proc2')
        def work2 = new WorkflowDef(name:'work2')
        meta2.addDefinition(proc2, func2, work2)

        // defs in the third script imported in a separate namespace
        def func3 = new FunctionDef(name: 'func3')
        def proc3 = new ProcessDef(script2, Mock(Closure), 'proc3')
        def work3 = new WorkflowDef(name:'work3')
        meta3.addDefinition(proc3, func3, work3)

        when:
        meta1.addModule(meta2, null, null)
        meta1.addModule(meta3, 'proc3', 'my_process')
        meta1.addModule(meta3, 'work3', null)

        then:
        meta1.getDefinitions() as Set == [proc1, func1, work1] as Set
        meta1.getComponent('proc1') == proc1
        meta1.getComponent('func1') == func1
        meta1.getComponent('work1') == work1

        then:
        // from root namespace
        meta1.getComponent('proc2') == proc2
        meta1.getComponent('func2') == func2
        meta1.getComponent('work2') == work2

        then:
        // other namespace are not reachable
        meta1.getComponent('xxx') == null
        meta1.getComponent('func3') == null
        meta1.getComponent('proc3') == null
        meta1.getComponent('work3') == work3
        meta1.getComponent('my_process') instanceof ProcessDef
        meta1.getComponent('my_process').name == 'my_process'

//        then:
//        meta1.getProcessNames() == ['proc1','proc2','my_process'] as Set
    }


    def 'should add a module component' () {
        given:
        def meta = Spy(ScriptMeta)
        def comp1 = Mock(ComponentDef)
        def comp2 = Mock(ComponentDef)

        // should a component to imports
        when:
        meta.addModule0(comp1)
        then:
        2 * comp1.getName() >> 'foo'
        1 * meta.getComponent('foo') >> null
        meta.@imports.get('foo') == comp1

        // should a component to imports with alias name
        when:
        meta.@imports.clear()
        meta.addModule0(comp1, 'bar')
        then:
        1 * comp1.getName() >> 'foo'
        1 * meta.getComponent('bar') >> null
        1 * comp1.cloneWithName('bar') >> comp2
        meta.@imports.get('bar') == comp2

    }

    def 'should throw a duplicate process name exception' () {
        given:
        def meta = Spy(ScriptMeta)
        def comp1 = Mock(ComponentDef)

        when:
        meta.@imports.clear()
        meta.addModule0(comp1)
        then:
        1 * comp1.getName() >> 'foo'
        1 * meta.getComponent('foo') >> comp1

        thrown(DuplicateModuleIncludeException)
    }
}
