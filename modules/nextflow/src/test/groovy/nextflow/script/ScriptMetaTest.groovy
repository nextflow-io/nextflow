/*
 * Copyright 2013-2024, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.script

import java.nio.file.Files

import groovy.transform.InheritConstructors
import nextflow.NF
import nextflow.exception.DuplicateModuleFunctionException
import test.Dsl2Spec
import test.TestHelper

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

    def setupSpec(){
        NF.init()
    }

    def 'should return all defined names' () {

        given:
        def script = new FooScript(new ScriptBinding())

        def proc1 = new ProcessDef(script, Mock(Closure), 'proc1')
        def proc2 = new ProcessDef(script, Mock(Closure), 'proc2')
        def func1 = new FunctionDef(name: 'func1', alias: 'func1')
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
        def func1 = new FunctionDef(name: 'func1', alias: 'func1')
        def proc1 = new ProcessDef(script1, Mock(Closure), 'proc1')
        def work1 = new WorkflowDef(name:'work1')
        meta1.addDefinition(proc1, func1, work1)

        // defs in the second script imported in the root namespace
        def func2 = new FunctionDef(name: 'func2', alias: 'func2')
        def proc2 = new ProcessDef(script2, Mock(Closure), 'proc2')
        def work2 = new WorkflowDef(name:'work2')
        meta2.addDefinition(proc2, func2, work2)

        // defs in the third script imported in a separate namespace
        def func3 = new FunctionDef(name: 'func3', alias: 'func3')
        def proc3 = new ProcessDef(script2, Mock(Closure), 'proc3')
        def work3 = new WorkflowDef(name:'work3')
        meta3.addDefinition(proc3, func3, work3)

        when:
        meta1.addModule(meta2, 'func2', null)
        meta1.addModule(meta2, 'proc2', null)
        meta1.addModule(meta2, 'work2', null)
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
        and:
        meta1.getComponent('work3') == work3
        meta1.getComponent('my_process') instanceof ProcessDef
        meta1.getComponent('my_process').name == 'my_process'

        then:
        meta1.getProcessNames() == ['proc1','proc2','my_process'] as Set
        and:
        meta1.getWorkflowNames() == ['work1', 'work2', 'work3'] as Set
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
        meta.@imports.get('foo') == comp1

        // should a component to imports with alias name
        when:
        meta.@imports.clear()
        meta.addModule0(comp1, 'bar')
        then:
        1 * comp1.getName() >> 'foo'
        1 * comp1.cloneWithName('bar') >> comp2
        meta.@imports.get('bar') == comp2

    }

    def 'should not throw a duplicate process name exception' () {
        given:
        def meta = Spy(ScriptMeta)
        def comp1 = Mock(ComponentDef)

        when:
        meta.@imports.clear()
        meta.addModule0(comp1)
        then:
        2 * comp1.getName() >> 'foo'
    }

    def 'should get module bundle' () {
        given:
        def folder = TestHelper.createInMemTempDir()
        and:
        def mod = folder.resolve('mod1'); mod.mkdir()
        mod.resolve('resources').mkdir()
        and:
        def scriptPath = mod.resolve('main.nf')
        def dockerPath = mod.resolve('Dockerfile')
        Files.createFile(scriptPath)
        Files.createFile(dockerPath)
        Files.createFile(mod.resolve('resources/foo.txt'))
        Files.createFile(mod.resolve('resources/bar.txt'))
        and:
        def meta = new ScriptMeta(scriptPath: scriptPath)

        when:
        def bundle = meta.getModuleBundle()
        then:
        bundle.dockerfile == dockerPath
        bundle.getEntries() == ['foo.txt', 'bar.txt'] as Set

    }

    def 'should throw duplicate name exception' () {

        given:
        def script1 = new FooScript(new ScriptBinding())
        def script2 = new FooScript(new ScriptBinding())
        def meta1 = new ScriptMeta(script1)
        def meta2 = new ScriptMeta(script2)

        // import module into main script
        def func2 = new FunctionDef(name: 'func1', alias: 'func1')
        def proc2 = new ProcessDef(script2, Mock(Closure), 'proc1')
        def work2 = new WorkflowDef(name: 'work1')
        meta2.addDefinition(proc2, func2, work2)

        meta1.addModule(meta2, 'func1', null)
        meta1.addModule(meta2, 'proc1', null)
        meta1.addModule(meta2, 'work1', null)

        // attempt to define duplicate components in main script
        def func1 = new FunctionDef(name: 'func1', alias: 'func1')
        def proc1 = new ProcessDef(script1, Mock(Closure), 'proc1')
        def work1 = new WorkflowDef(name: 'work1')

        when:
        meta1.addDefinition(func1)
        then:
        def e = thrown(DuplicateModuleFunctionException)
        e.message.contains "A function named 'func1' is already defined"

        when:
        meta1.addDefinition(proc1)
        then:
        e = thrown(DuplicateModuleFunctionException)
        e.message.contains "A process named 'proc1' is already defined"

        when:
        meta1.addDefinition(work1)
        then:
        e = thrown(DuplicateModuleFunctionException)
        e.message.contains "A workflow named 'work1' is already defined"
    }
}
