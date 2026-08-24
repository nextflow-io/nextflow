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

import java.lang.management.ManagementFactory

import com.sun.management.OperatingSystemMXBean
import nextflow.Session
import nextflow.SysEnv
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.ExecutorConfig
import nextflow.util.MemoryUnit
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LocalPollingMonitorTest extends Specification {

    def 'should allocated and free resources' () {

        given:
        def _20_GB = MemoryUnit.of('20GB').toBytes()
        def session = Mock(Session)
        def monitor = new LocalPollingMonitor(
                cpus: 10,
                capacity: 20,
                memory: _20_GB,
                session: session,
                name: 'local',
                pollInterval: 100
        )

        def task = new TaskRun()
        task.config = new TaskConfig(cpus: 3, memory: MemoryUnit.of('2GB'))
        def handler = Mock(TaskHandler)
        handler.getTask() >> { task }

        expect:
        monitor.availCpus == 10
        monitor.capacity == 20
        monitor.availMemory == _20_GB
        monitor.maxCpus == 10
        monitor.maxMemory == _20_GB

        when:
        monitor.submit(handler)
        then:
        session.notifyTaskSubmit(handler) >> null
        monitor.getRunningQueue().size()==1
        monitor.availCpus == 7
        monitor.availMemory == MemoryUnit.of('18GB').toBytes()
        monitor.maxCpus == 10
        monitor.maxMemory == _20_GB

        when:
        monitor.remove(handler)
        then:
        monitor.getRunningQueue().size()==0
        monitor.availCpus == 10
        monitor.availMemory == _20_GB
        monitor.maxCpus == 10
        monitor.maxMemory == _20_GB

    }


    def 'validate canSubmit method' () {

        given:
        def _20_GB = MemoryUnit.of('20GB').toBytes()
        def session = Mock(Session)
        def monitor = new LocalPollingMonitor(
                cpus: 10,
                capacity: 10,
                memory: _20_GB,
                session: session,
                name: 'local',
                pollInterval: 100
        )

        def task = new TaskRun()
        task.config = new TaskConfig(cpus: 4, memory: MemoryUnit.of('8GB'))
        def handler = Mock(TaskHandler)
        handler.getTask() >> { task }
        handler.canForkProcess() >> true
        handler.isReady() >> true

        expect:
        monitor.canSubmit(handler) == true

        when:
        monitor.submit(handler)
        then:
        1 * handler.submit() >> null
        1 * session.notifyTaskSubmit(handler) >> null
        and:
        monitor.canSubmit(handler) == true
        monitor.availCpus == 6
        monitor.availMemory == MemoryUnit.of('12GB').toBytes()

        when:
        monitor.submit(handler)
        then:
        1 * handler.submit() >> null
        1 * session.notifyTaskSubmit(handler) >> null
        and:
        monitor.canSubmit(handler) == false
        monitor.availCpus == 2
        monitor.availMemory == MemoryUnit.of('4GB').toBytes()

    }

    def 'validate canSubmit method for one cpus' () {

        given:
        def _20_GB = MemoryUnit.of('20GB').toBytes()
        def session = Mock(Session)
        def monitor = new LocalPollingMonitor(
                cpus: 1,
                capacity: 1,
                memory: _20_GB,
                session: session,
                name: 'local',
                pollInterval: 100
        )

        def task = new TaskRun()
        task.config = new TaskConfig(cpus: 1, memory: MemoryUnit.of('8GB'))
        def handler = Mock(TaskHandler)
        handler.getTask() >> { task }
        handler.canForkProcess() >> true
        handler.isReady() >> true

        expect:
        monitor.canSubmit(handler) == true
        monitor.availCpus == 1

        when:
        monitor.submit(handler)
        then:
        1 * handler.submit() >> null
        1 * session.notifyTaskSubmit(handler) >> null
        and:
        monitor.canSubmit(handler) == false
        monitor.availCpus == 0
    }

    def 'should throw an exception for missing cpus' () {

        given:
        def _20_GB = MemoryUnit.of('20GB').toBytes()
        def session = new Session()
        def monitor = new LocalPollingMonitor(
                cpus: 10,
                capacity: 20,
                memory: _20_GB,
                session: session,
                name: 'local',
                pollInterval: 100
        )

        def task = new TaskRun()
        task.config = new TaskConfig(cpus: 12)
        def handler = Mock(TaskHandler)
        handler.getTask() >> { task }

        when:
        monitor.canSubmit(handler)
        then:
        def e1 = thrown(ProcessUnrecoverableException)
        e1.message == 'Process requirement exceeds available CPUs -- req: 12; avail: 10'


    }

    def 'should throw an exception for missing memory' () {

        given:
        def _20_GB = MemoryUnit.of('20GB').toBytes()
        def session = new Session()
        def monitor = new LocalPollingMonitor(
                cpus: 10,
                capacity: 20,
                memory: _20_GB,
                session: session,
                name: 'local',
                pollInterval: 100
        )

        def task = new TaskRun()
        task.config = new TaskConfig(memory: MemoryUnit.of('22GB'))
        def handler = Mock(TaskHandler)
        handler.getTask() >> { task }

        when:
        monitor.canSubmit(handler)
        then:
        def e2 = thrown(ProcessUnrecoverableException)
        e2.message == 'Process requirement exceeds available memory -- req: 22 GB; avail: 20 GB'

    }

    def 'should throw an exception for missing accelerators' () {

        given:
        SysEnv.push(CUDA_VISIBLE_DEVICES: '0,1,2,3')
        and:
        def _20_GB = MemoryUnit.of('20GB').toBytes()
        def session = new Session()
        def monitor = new LocalPollingMonitor(
                cpus: 10,
                capacity: 20,
                memory: _20_GB,
                session: session,
                name: 'local',
                pollInterval: 100
        )
        and:
        def task = new TaskRun()
        task.config = new TaskConfig(accelerator: 8)
        def handler = Mock(TaskHandler)
        handler.getTask() >> { task }

        when:
        monitor.canSubmit(handler)
        then:
        def e2 = thrown(ProcessUnrecoverableException)
        e2.message == 'Process requirement exceeds available accelerators -- req: 8; avail: 4'

        cleanup:
        SysEnv.pop()

    }

    def 'should get the number of cpus' () {

        given:
        def OS = (OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean()

        when:
        def config1 = new ExecutorConfig([:])
        then:
        LocalPollingMonitor.configCpus(config1,'local') == OS.getAvailableProcessors()

        when:
        def config2 = new ExecutorConfig([cpus: 100])
        then:
        LocalPollingMonitor.configCpus(config2,'local') == 100

        when:
        def config3 = new ExecutorConfig(['$local': [cpus: 100]])
        then:
        LocalPollingMonitor.configCpus(config3,'local') == 100

        when:
        def config4 = new ExecutorConfig(['$sge': [cpus: 100]])
        then:
        LocalPollingMonitor.configCpus(config4,'local') == OS.getAvailableProcessors()

    }

    def 'should get the amount of mem' () {

        given:
        def OS = (OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean()
        def _10_GB = MemoryUnit.of('10 GB').toBytes()

        when:
        def config1 = new ExecutorConfig([:])
        then:
        LocalPollingMonitor.configMem(config1,'local') == OS.getTotalPhysicalMemorySize()

        when:
        def config2 = new ExecutorConfig([memory: '10 GB'])
        then:
        LocalPollingMonitor.configMem(config2,'local') == _10_GB

        when:
        def config3 = new ExecutorConfig(['$local': [memory: _10_GB]])
        then:
        LocalPollingMonitor.configMem(config3,'local') == _10_GB

        when:
        def config4 = new ExecutorConfig(['$sge': [memory: '1 GB']])
        then:
        LocalPollingMonitor.configMem(config4,'local') == OS.getTotalPhysicalMemorySize()

    }
}
