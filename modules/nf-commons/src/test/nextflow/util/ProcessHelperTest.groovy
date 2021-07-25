package nextflow.util

import java.lang.management.ManagementFactory

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessHelperTest extends Specification {

    def 'should get process pid' () {

        given:
        def process = new ProcessBuilder().command(['bash','-c','echo $$']).start()

        when:
        def pid = ProcessHelper.pid(process)
        and:
        process.waitFor()
        then:
        process.text.trim() == pid.toString()
    }

    def 'should return self pid' () {
        when:
        def pid = ProcessHelper.selfPid()
        then:
        pid > 0
        pid == Long.parseLong(ManagementFactory.getRuntimeMXBean().getName().split("@")[0])
    }
}
