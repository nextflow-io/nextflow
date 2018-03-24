/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.processor
import nextflow.Session
import nextflow.executor.AwsBatchExecutor
import nextflow.executor.CondorExecutor
import nextflow.executor.Executor
import nextflow.executor.LocalExecutor
import nextflow.executor.LsfExecutor
import nextflow.executor.NopeExecutor
import nextflow.executor.PbsExecutor
import nextflow.executor.SgeExecutor
import nextflow.executor.SlurmExecutor
import nextflow.k8s.K8sExecutor
import nextflow.script.BaseScript
import nextflow.script.ScriptType
import nextflow.util.ServiceName
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessFactoryTest extends Specification {


    def 'test loadExecutor' () {

        setup:
        def factory = [:] as ProcessFactory
        factory.executorsMap [ 'x' ] = XExecutor.name   // // <-- this is loaded by the name
        expect:
        factory.loadExecutorClass(null) == LocalExecutor
        factory.loadExecutorClass('local') == LocalExecutor
        factory.loadExecutorClass('sge') == SgeExecutor
        factory.loadExecutorClass('oge') == SgeExecutor
        factory.loadExecutorClass('uge') == SgeExecutor
        factory.loadExecutorClass('lsf') == LsfExecutor
        factory.loadExecutorClass('pbs') == PbsExecutor
        factory.loadExecutorClass('slurm') == SlurmExecutor
        factory.loadExecutorClass('condor') == CondorExecutor
        factory.loadExecutorClass('k8s') == K8sExecutor
        factory.loadExecutorClass('awsbatch') == AwsBatchExecutor
        factory.loadExecutorClass('AwsBatch') == AwsBatchExecutor
        factory.loadExecutorClass('aws-batch') == AwsBatchExecutor
        factory.loadExecutorClass('x') == XExecutor  // <-- this is loaded by the name

        when:
        factory.loadExecutorClass('xyz')
        then:
        thrown(IllegalArgumentException)

    }

    def testLoader() {
        when:
        def factory = new ProcessFactory(Mock(BaseScript), Mock(Session) )
        then:
        factory.loadExecutorClass('sge') == SgeExecutor
    }


    def testSupportType() {

        setup:
        def factory = [:] as ProcessFactory

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


    def testGetExecName() {
        expect:
        ProcessFactory.findNameByClass( NopeExecutor ) == 'nope'
        ProcessFactory.findNameByClass( SgeExecutor ) == 'sge'
        ProcessFactory.findNameByClass(XExecutor) == 'my_fancy_name'
    }


    def 'should apply config setting for a process label' () {
        given:
        def factory = new ProcessFactory()
        def settings = [
                'withLabel:short'  : [ cpus: 1, time: '1h'],
                'withLabel:!short' : [ cpus: 32, queue: 'cn-long'],
                'withLabel:foo'    : [ cpus: 2 ],
                'withLabel:foo|bar': [ disk: '100GB' ],
                'withLabel:gpu.+'  : [ cpus: 4 ],
        ]

        when:
        def process = new ProcessConfig([:])
        factory.applyConfigForLabel('any', process, 'withLabel:', 'short', settings)
        then:
        process.cpus == 1
        process.time == '1h'
        process.size() == 2

        when:
        process = new ProcessConfig([:])
        factory.applyConfigForLabel('any', process, 'withLabel:', 'long', settings)
        then:
        process.cpus == 32
        process.queue == 'cn-long'
        process.size() == 2

        when:
        process = new ProcessConfig([:])
        factory.applyConfigForLabel('any', process, 'withLabel:', 'foo', settings)
        then:
        process.cpus == 2
        process.disk == '100GB'
        process.queue == 'cn-long'
        process.size() == 3

        when:
        process = new ProcessConfig([:])
        factory.applyConfigForLabel('any', process, 'withLabel:', 'bar', settings)
        then:
        process.cpus == 32
        process.disk == '100GB'
        process.queue == 'cn-long'
        process.size() == 3

        when:
        process = new ProcessConfig([:])
        factory.applyConfigForLabel('any', process, 'withLabel:', 'gpu-1', settings)
        then:
        process.cpus == 4
        process.queue == 'cn-long'
        process.size() == 2
        
    }


    def 'should apply config setting for a process name' () {
        given:
        def factory = new ProcessFactory()
        def settings = [
                'withName:alpha'        : [ cpus: 1, time: '1h'],
                'withName:delta'        : [ cpus: 2 ],
                'withName:delta|gamma'  : [ disk: '100GB' ],
                'withName:omega.+'      : [ cpus: 4 ],
        ]

        when:
        def process = new ProcessConfig([:])
        factory.applyConfigForLabel('xx', process, 'withName:', 'any', settings)
        then:
        process.size() == 0

        when:
        process = new ProcessConfig([:])
        factory.applyConfigForLabel('alpha', process, 'withName:', 'any', settings)
        then:
        process.cpus == 1
        process.time == '1h'
        process.size() == 2

        when:
        process = new ProcessConfig([:])
        factory.applyConfigForLabel('delta', process, 'withName:', 'any', settings)
        then:
        process.cpus == 2
        process.disk == '100GB'
        process.size() == 2

        when:
        process = new ProcessConfig([:])
        factory.applyConfigForLabel('gamma', process, 'withName:', 'any', settings)
        then:
        process.disk == '100GB'
        process.size() == 1

        when:
        process = new ProcessConfig([:])
        factory.applyConfigForLabel('omega_x', process,'withName:', 'any', settings)
        then:
        process.cpus == 4
        process.size() == 1
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