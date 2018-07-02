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

import java.lang.management.ManagementFactory

import com.sun.management.OperatingSystemMXBean
import nextflow.Session
import nextflow.exception.ProcessUnrecoverableException
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
        e1.message == 'Process requirement exceed available CPUs -- req: 12; avail: 10'


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
        e2.message == 'Process requirement exceed available memory -- req: 22 GB; avail: 20 GB'

    }

    def 'should get the number of cpus' () {

        given:
        def OS = (OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean()

        when:
        def session1 = new Session()
        then:
        LocalPollingMonitor.configCpus(session1,'local') == OS.getAvailableProcessors()

        when:
        def session2 = new Session([executor: [cpus: 100]])
        then:
        LocalPollingMonitor.configCpus(session2,'local') == 100

        when:
        def session3 = new Session([executor: ['$local': [cpus: 100]]])
        then:
        LocalPollingMonitor.configCpus(session3,'local') == 100

        when:
        def session4 = new Session([executor: ['$sge': [cpus: 100]]])
        then:
        LocalPollingMonitor.configCpus(session4,'local') == OS.getAvailableProcessors()

    }

    def 'should get the amount of mem' () {

        given:
        def OS = (OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean()
        def _10_GB = MemoryUnit.of('10 GB').toBytes()

        when:
        def session1 = new Session()
        then:
        LocalPollingMonitor.configMem(session1,'local') == OS.getTotalPhysicalMemorySize()

        when:
        def session2 = new Session([executor: [memory: '10 GB']])
        then:
        LocalPollingMonitor.configMem(session2,'local') == _10_GB

        when:
        def session3 = new Session([executor: ['$local': [memory: _10_GB]]])
        then:
        LocalPollingMonitor.configMem(session3,'local') == _10_GB

        when:
        def session4 = new Session([executor: ['$sge': [memory: '1 GB']]])
        then:
        LocalPollingMonitor.configMem(session4,'local') == OS.getTotalPhysicalMemorySize()

    }
}
