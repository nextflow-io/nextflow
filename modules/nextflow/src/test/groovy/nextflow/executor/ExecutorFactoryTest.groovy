/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.executor

import spock.lang.Specification

import nextflow.Session
import nextflow.k8s.K8sExecutor
import nextflow.script.ProcessConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.script.ScriptType
import nextflow.script.BodyDef
import nextflow.util.ServiceName
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ExecutorFactoryTest extends Specification {


    def 'should load executor class' () {

        setup:
        def factory = new ExecutorFactory()
        factory.executorsMap [ 'x' ] = XExecutor.name   // // <-- this is loaded by the name
        expect:
        factory.getExecutorClass(null) == LocalExecutor
        factory.getExecutorClass('local') == LocalExecutor
        factory.getExecutorClass('sge') == SgeExecutor
        factory.getExecutorClass('oge') == SgeExecutor
        factory.getExecutorClass('uge') == SgeExecutor
        factory.getExecutorClass('lsf') == LsfExecutor
        factory.getExecutorClass('pbs') == PbsExecutor
        factory.getExecutorClass('slurm') == SlurmExecutor
        factory.getExecutorClass('condor') == CondorExecutor
        factory.getExecutorClass('k8s') == K8sExecutor
        factory.getExecutorClass('x') == XExecutor  // <-- this is loaded by the name

        when:
        factory.getExecutorClass('xyz')
        then:
        thrown(IllegalArgumentException)

    }

    def 'should not duplicate executor instance' () {
        given:
        def NAME = "SLURM"
        def config = Mock(ProcessConfig)
        def session = Mock(Session)
        def script = Mock(BodyDef)
        script.type >> ScriptType.SCRIPTLET
        def factory = Spy(ExecutorFactory)
        def clazz = SlurmExecutor.class
        def executor = Mock(SlurmExecutor)

        when:
        def result = factory.getExecutor('proc_a', config, script, session)
        then:
        1 * factory.getExecutorName(config, session) >> NAME
        1 * factory.getExecutorClass(NAME) >> clazz
        1 * factory.isTypeSupported(ScriptType.SCRIPTLET, clazz) >> true
        1 * factory.createExecutor(clazz, NAME, session) >> executor
        result.is(executor)

        when:
        factory.executors.put(clazz, executor)
        result = factory.getExecutor('proc_b', config, script, session)
        then:
        1 * factory.getExecutorName(config, session) >> NAME
        1 * factory.getExecutorClass(NAME) >> clazz
        1 * factory.isTypeSupported(ScriptType.SCRIPTLET, clazz) >> true
        0 * factory.createExecutor(clazz, NAME, session)
        result.is(executor)

    }
    

    def 'should check type supported'() {

        setup:
        def factory = new ExecutorFactory()

        when:
        factory.isTypeSupported(ScriptType.GROOVY, 'xxx')
        then:
        thrown(IllegalArgumentException)

        expect:
        // by default only script are supported
        factory.isTypeSupported(ScriptType.SCRIPTLET, Executor,)
        !factory.isTypeSupported(ScriptType.GROOVY, Executor)
        // same for grid
        factory.isTypeSupported(ScriptType.SCRIPTLET, SgeExecutor)
        !factory.isTypeSupported(ScriptType.GROOVY, SgeExecutor)

        // local supports both                               x
        factory.isTypeSupported(ScriptType.SCRIPTLET, LocalExecutor)
        factory.isTypeSupported(ScriptType.GROOVY, LocalExecutor)

        // repeat for instances
        factory.isTypeSupported(ScriptType.SCRIPTLET, new SgeExecutor() )
        !factory.isTypeSupported(ScriptType.GROOVY, new SgeExecutor())
        factory.isTypeSupported(ScriptType.SCRIPTLET, new LocalExecutor())
        factory.isTypeSupported(ScriptType.GROOVY, new LocalExecutor())
    }


    def 'should return executor name'() {
        expect:
        ExecutorFactory.findNameByClass( NopeExecutor ) == 'nope'
        ExecutorFactory.findNameByClass( SgeExecutor ) == 'sge'
        ExecutorFactory.findNameByClass(XExecutor) == 'my_fancy_name'
    }
}


@ServiceName('my_fancy_name')
class XExecutor extends Executor {
    @Override
    protected TaskMonitor createTaskMonitor() {
        return null
    }

    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        return null
    }
}
