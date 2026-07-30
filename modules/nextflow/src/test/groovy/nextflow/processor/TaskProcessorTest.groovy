/*
 * Copyright 2013-2026, Seqera Labs
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
import java.nio.file.Path
import java.nio.file.Paths

import groovyx.gpars.agent.Agent
import nextflow.Session
import nextflow.executor.Executor
import nextflow.executor.NopeExecutor
import nextflow.file.FileHolder
import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import nextflow.script.ScriptMeta
import nextflow.script.ScriptType
import nextflow.script.bundle.ResourcesBundle
import nextflow.util.CacheHelper
import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskProcessorTest extends Specification {

    def createProcessor(String name, Session session) {
        return new DummyProcessor(name, session, Mock(BaseScript), new ProcessConfig([:]))
    }

    static class DummyProcessor extends TaskProcessor {

        DummyProcessor(String name, Session session, BaseScript script, ProcessConfig config) {
            super(name, new NopeExecutor(session: session), session, script, config, new BodyDef({}, '..'))
        }

        @Override
        protected void createOperator() { }
    }


    def 'should evaluate environment variables'() {

        setup:
        def home = Files.createTempDirectory('test')
        def binFolder = home.resolve('bin')
        binFolder.mkdirs()

        when:
        def session = new Session([env: [X:"1", Y:"2", Z:null, W:'']])
        session.setBaseDir(home)
        def processor = createProcessor('task1', session)
        def builder = new ProcessBuilder()
        builder.environment().putAll( processor.getProcessEnvironment() )
        then:
        noExceptionThrown()
        builder.environment().X == '1'
        builder.environment().Y == '2'
        builder.environment().PATH == "\$PATH:${binFolder.toString()}"
        !builder.environment().containsKey('Z')
        // explicit empty-string values are retained (unlike null values) -- see #5722
        builder.environment().W == ''

        when:
        session = new Session([env: [X:"1", Y:"2", PATH:'/some']])
        session.setBaseDir(home)
        processor = createProcessor('task1', session)
        builder = new ProcessBuilder()
        builder.environment().putAll( processor.getProcessEnvironment() )
        then:
        noExceptionThrown()
        builder.environment().X == '1'
        builder.environment().Y == '2'
        builder.environment().PATH == "/some:${binFolder.toString()}"

        cleanup:
        home.deleteDir()

    }

    @Unroll
    def 'should resolve module bundle for entry script when running as module #desc' () {
        given:
        def folder = Files.createTempDirectory('test')
        def mod = folder.resolve('mod1'); mod.mkdir()
        def bin = mod.resolve('resources/usr/bin'); bin.mkdirs()
        def scriptPath = mod.resolve('main.nf'); Files.createFile(scriptPath)
        Files.createFile(bin.resolve('echo.sh'))
        and:
        def script = Mock(BaseScript)
        def meta = Mock(ScriptMeta) {
            getScriptPath() >> scriptPath
            isModule() >> IS_MODULE
            getModuleBundle() >> ResourcesBundle.scan(mod.resolve('resources'))
        }
        and:
        def session = Mock(Session) {
            getConfig() >> [:]
            isModuleRun() >> IS_MODULE_RUN
        }
        def executor = Mock(Executor) {}
        def processor = Spy(TaskProcessor, constructorArgs: [[session:session, executor:executor]])
        processor.getOwnerScript() >> script

        when:
        ResourcesBundle bundle
        GroovyMock(ScriptMeta, global: true)
        ScriptMeta.get(script) >> meta
        bundle = processor.getModuleBundle()

        then:
        (bundle != null) == EXPECTED

        cleanup:
        folder?.deleteDir()

        where:
        desc                                | IS_MODULE | IS_MODULE_RUN | EXPECTED
        '(included module)'                 | true      | false         | true
        '(nextflow module run entry)'       | false     | true          | true
        '(plain entry, no module run flag)' | false     | false         | false
    }

    @Unroll
    def 'should add module bin paths to task env' () {
        given:
        def session = Mock(Session) { getConfig() >> [:] }
        def executor = Mock(Executor) { getBinDir() >> Path.of('/project/bin')}
        and:
        TaskProcessor processor = Spy(TaskProcessor, constructorArgs: [[session:session, executor:executor]])
        and:
        when:
        def result = processor.getProcessEnvironment()
        then:
        session.enableModuleBinaries() >> MODULE_BIN
        processor.getModuleBundle() >> Mock(ResourcesBundle)  { getBinDirs() >> [Path.of('/foo'), Path.of('/bar')] }
        processor.isLocalWorkDir() >> LOCAL
        and:
        result == EXPECTED

        where:
        LOCAL   | MODULE_BIN    | EXPECTED
        false   | false         | [:]
        true    | false         | [PATH:'$PATH:/project/bin']
        true    | true          | [PATH:'$PATH:/foo:/bar:/project/bin']
    }

    def 'should fetch interpreter from shebang line'() {

        when:
        def script =
            '''
            #!/bin/perl
            do this
            do that
            '''
        def i = TaskProcessor.fetchInterpreter(script.stripIndent().trim())
        then:
        i == '/bin/perl'

        when:
        i = TaskProcessor.fetchInterpreter('do this')
        then:
        i == null
    }

    def 'should stage path'() {

        when:
        def p1 = Paths.get('/home/data/source.file')
        def path = FileHolder.get(p1, 'target.file')

        then:
        path.storePath == p1
        path.stageName == 'target.file'

    }


    def 'should return task hash log'() {

        when:
        def h = CacheHelper.hasher('x').hash()
        def task = new TaskRun(hash:h)
        then:
        task.getHashLog() == '76/9f897d'

    }

    def 'should update agent state'() {
        when:
        def state = new Agent<StateObj>(new StateObj())
        int i = 0
        state.addListener { a, b -> i++ }

        state.update { StateObj it ->  it.incSubmitted()  }
        state.update { StateObj it ->  it.incCompleted() }
        state.update  { StateObj it ->  it.poison()  }
        state.await()
        then:
        state.val.finished
        i == 3

        when:
        state = new Agent<StateObj>(new StateObj())
        state.update { StateObj it ->  it.incSubmitted()  }
        state.update  { StateObj it ->  it.poison()  }
        state.await()
        then:
        !state.val.finished

    }

    def 'should update agent state 2'() {

        when:
        def agent = new Agent<List>([])
        int i = 0
        agent.addListener { a, b -> println ">>: $a -- $b"; i++ }

        agent << { it.add(1); (this as Agent).updateValue(it) }
        agent << { it.add(2); (this as Agent).updateValue(it) }
        agent.await()
        then:
        agent.val == [1,2]

    }

    def 'should export env vars' () {

        given:
        def env

        when:
        env = TaskProcessor.bashEnvironmentScript([FOO:'hola',BAR:'ciao mondo'])
        then:
        env ==  '''
                export FOO="hola"
                export BAR="ciao mondo"
                '''
                .stripIndent().leftTrim()

        when:
        env = TaskProcessor.bashEnvironmentScript([PATH: 'foo:$PATH', HOLA: 'one|two'], true)
        then:
        env == '''\
            export PATH="foo:\\$PATH"
            export HOLA="one|two"
            '''.stripIndent()
        env.charAt(env.size()-1) == '\n' as char

        when:
        env = TaskProcessor.bashEnvironmentScript([FOO:null, BAR:''])
        then:
        env == "export FOO=''\nexport BAR=''\n"

    }

    def 'should bind fair outputs' () {
        given:
        def processor = Spy(TaskProcessor)
        processor.@config = Mock(ProcessConfig)
        processor.@isFair0 = true
        and:
        def task3 = Mock(TaskRun) { getIndex()>>3 }
        and:
        def task2 = Mock(TaskRun) { getIndex()>>2 }
        and:
        def task1 = Mock(TaskRun) { getIndex()>>1 }
        and:
        def task5 = Mock(TaskRun) { getIndex()>>5 }
        and:
        def task4 = Mock(TaskRun) { getIndex()>>4 }

        when:
        processor.fairBindOutputs0(task3)
        then:
        processor.@fairBuffers[2] == task3
        0 * processor.bindOutputs0(_)

        when:
        processor.fairBindOutputs0(task2)
        then:
        processor.@fairBuffers[1] == task2
        0 * processor.bindOutputs0(_)

        when:
        processor.fairBindOutputs0(task5)
        then:
        processor.@fairBuffers[4] == task5
        0 * processor.bindOutputs0(_)

        when:
        processor.fairBindOutputs0(task1)
        then:
        1 * processor.bindOutputs0(task1)
        then:
        1 * processor.bindOutputs0(task2)
        then:
        1 * processor.bindOutputs0(task3)
        and:
        processor.@fairBuffers.size() == 2
        processor.@fairBuffers[0] == null
        processor.@fairBuffers[1] == task5

        when:
        processor.fairBindOutputs0(task4)
        then:
        1 * processor.bindOutputs0(task4)
        then:
        1 * processor.bindOutputs0(task5)
        then:
        processor.@fairBuffers.size()==0
    }

    def 'should create a task preview' () {
        given:
        def config = new ProcessConfig([cpus: 10, memory: '100 GB'])
        def EXEC = Mock(Executor) { getName()>>'exec-name'}
        def BODY = Mock(BodyDef) { getType()>>ScriptType.SCRIPTLET }
        def processor = new TaskProcessor(config: config, name: 'proc-name', executor: EXEC, taskBody: BODY)

        when:
        def result = processor.createTaskPreview()
        then:
        result.config.process == 'proc-name'
        result.config.executor == 'exec-name'
        result.config.getCpus() == 10
        result.config.getMemory() == MemoryUnit.of('100 GB')
    }

}
