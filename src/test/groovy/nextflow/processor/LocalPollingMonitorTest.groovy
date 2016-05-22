package nextflow.processor

import java.lang.management.ManagementFactory

import com.sun.management.OperatingSystemMXBean
import nextflow.Session
import nextflow.exception.ProcessNotRecoverableException
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
        def session = new Session()
        def monitor = new LocalPollingMonitor(
                cpus: 10,
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
        monitor.availMemory == _20_GB
        monitor.maxCpus == 10
        monitor.maxMemory == _20_GB

        when:
        monitor.submit(handler)
        then:
        monitor.getPollingQueue().size()==1
        monitor.availCpus == 7
        monitor.availMemory == MemoryUnit.of('18GB').toBytes()
        monitor.maxCpus == 10
        monitor.maxMemory == _20_GB

        when:
        monitor.remove(handler)
        then:
        monitor.getPollingQueue().size()==0
        monitor.availCpus == 10
        monitor.availMemory == _20_GB
        monitor.maxCpus == 10
        monitor.maxMemory == _20_GB

    }


    def 'validate canSubmit method' () {

        given:
        def _20_GB = MemoryUnit.of('20GB').toBytes()
        def session = new Session()
        def monitor = new LocalPollingMonitor(
                cpus: 10,
                memory: _20_GB,
                session: session,
                name: 'local',
                pollInterval: 100
        )

        def task = new TaskRun()
        task.config = new TaskConfig(cpus: 4, memory: MemoryUnit.of('8GB'))
        def handler = Mock(TaskHandler)
        handler.getTask() >> { task }

        when:
        true
        then:
        monitor.canSubmit(handler)

        when:
        monitor.submit(handler)
        then:
        monitor.canSubmit(handler)
        monitor.availCpus == 6
        monitor.availMemory == MemoryUnit.of('12GB').toBytes()

        when:
        monitor.submit(handler)
        then:
        !monitor.canSubmit(handler)
        monitor.availCpus == 2
        monitor.availMemory == MemoryUnit.of('4GB').toBytes()


    }

    def 'should throw an exception for missing cpus' () {

        given:
        def _20_GB = MemoryUnit.of('20GB').toBytes()
        def session = new Session()
        def monitor = new LocalPollingMonitor(
                cpus: 10,
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
        def e1 = thrown(ProcessNotRecoverableException)
        e1.message == 'Process requirement exceed available CPUs -- req: 12; avail: 10'


    }

    def 'should throw an exception for missing memory' () {

        given:
        def _20_GB = MemoryUnit.of('20GB').toBytes()
        def session = new Session()
        def monitor = new LocalPollingMonitor(
                cpus: 10,
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
        def e2 = thrown(ProcessNotRecoverableException)
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
