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

package nextflow.processor

import java.nio.file.Files
import java.nio.file.Paths

import groovy.runtime.metaclass.ExtensionProvider
import groovy.runtime.metaclass.NextflowDelegatingMetaClass
import groovy.transform.InheritConstructors
import nextflow.Global
import nextflow.NF
import nextflow.NextflowMeta
import nextflow.Session
import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import nextflow.script.ScriptBinding
import nextflow.script.ScriptMeta
import nextflow.util.BlankSeparatedList
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskContextTest extends Specification {

    def setupSpec() {
        NF.init()
    }
    
    def 'should save and read TaskContext object' () {

        setup:
        def taskConfig = new ProcessConfig([:])
        def processor = [:] as TaskProcessor
        processor.metaClass.getTaskConfig = { taskConfig }
        processor.metaClass.getTaskBody = { new BodyDef(null,'source') }
        def str = 'Hola'
        def map = new TaskContext(processor, [:])
        map.alpha = 1
        map.beta = "${str}.txt"
        map.delta = new Duration('1day')
        map.micro = new MemoryUnit('100KB')
        map.file = Paths.get('Hola.txt')
        map.list = new BlankSeparatedList( Paths.get('A'), Paths.get('B'), Paths.get('C') )
        map.holder = 'just a string'
        map.context = [uno: 1, due: 'str', tre: "$str"]

        when:
        def buffer = map.serialize()
        def result = TaskContext.deserialize(processor, buffer)

        then:
        result.size() == 8
        result.alpha == 1
        result.beta == "${str}.txt"
        result.delta == new Duration('1day')
        result.micro == new MemoryUnit('100KB')
        result.file.equals( Paths.get('Hola.txt') )
        result.context == [uno: 1, due: 'str', tre: "$str"]
        result.list == new BlankSeparatedList( Paths.get('A'), Paths.get('B'), Paths.get('C') )
        result.holder == 'just a string'
        result.get('holder') == 'just a string'
        result.getHolder() instanceof Map
        result.getHolder().get('alpha') == 1

    }

    def 'should dehydrate rehydrate'() {

        setup:
        def bind = new ScriptBinding([x:1, y:2])
        def script = Mock(BaseScript) { getBinding() >> bind }
        and:
        def local = [p:3, q:4, path: Paths.get('some/path')]
        def delegate = new TaskContext( script, local, 'hola' )

        when:
        def bytes = delegate.dehydrate()
        def copy = TaskContext.rehydrate(bytes)

        then:
        delegate == copy
        delegate.getHolder() == copy.getHolder()
        copy.getHolder() == local

    }

    def 'should interpolate variables' () {

        given:
        GroovyShell shell

        when:
        shell = new GroovyShell(new Binding([foo: 'hello', params:[bar:'world']]))
        then:
        shell.evaluate('"$foo - ${params.bar}"').toString() == 'hello - world'

        when:
        def binding = new Binding(params:[bar:'world'])
        def local = [foo: 'hello']

        def script = Mock(Script)
        script.getBinding() >> binding

        def context = new TaskContext(script, local, 'hola')
        shell = new GroovyShell(new Binding(context))
        then:
        shell.evaluate('"$foo - ${params.bar}"').toString() == 'hello - world'

    }

    def 'should clone the context' () {
        given:
        def holder = [:]
        def context = new TaskContext(Mock(Script), holder, 'proc_1')

        when:
        context.xxx = 1
        then:
        // both the `context` and the `holder` have the same `xxx` property
        context.xxx == 1
        holder.xxx == 1

        when:
        def copy = context.clone()
        copy.yyy = 2
        then:
        // the `yyy` is not set the holder from where the `copy` context has been cloned
        copy.xxx == 1
        copy.yyy == 2
        holder.yyy == null

    }


    def 'should resolve absolute paths as template paths' () {
        given:
        def temp = Files.createTempDirectory('test')
        and:
        Global.session = Mock(Session) { getBaseDir() >> temp  }
        and:
        def holder = [:]
        def script = Mock(BaseScript)
        TaskContext context = Spy(TaskContext, constructorArgs: [script, holder, 'proc_1'])

        when:
        // an absolute path is specified
        def absolutePath = Paths.get('/some/template.txt')
        def result = context.template(absolutePath)
        then:
        // the path is returned
        result == absolutePath

        when:
        // an absolute string path is specified
        result = context.template('/some/template.txt')
        then:
        // the path is returned
        result == absolutePath

        when:
        // when a rel file path is matched against the project
        // template dir, it's returned as the template file
        temp.resolve('templates').mkdirs()
        temp.resolve('templates/foo.txt').text = 'echo hello'
        and:
        result = context.template('foo.txt')
        then:
        result == temp.resolve('templates/foo.txt')

        when:
        // it's a DSL2 module
        NextflowMeta.instance.enableDsl2()
        NextflowDelegatingMetaClass.provider = Mock(ExtensionProvider) { operatorNames() >> new HashSet<String>() }
        def meta = ScriptMeta.register(script)
        meta.setScriptPath(temp.resolve('modules/my-module/main.nf'))
        and:
        // the module has a nested `templates` directory
        temp.resolve('modules/my-module/templates').mkdirs()
        temp.resolve('modules/my-module/templates/foo.txt').text = 'echo Hola'
        and:
        // the template is a relative path
        result = context.template('foo.txt')
        then:
        // the path is resolved against the module templates
        result == temp.resolve('modules/my-module/templates/foo.txt')

        cleanup:
        NextflowDelegatingMetaClass.provider = null
        NextflowMeta.instance.disableDsl2()
        temp?.deleteDir()
    }


}

@InheritConstructors
class MockScript extends BaseScript {
    @Override Object runScript() { return null }
}
