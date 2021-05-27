package nextflow.util

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
}
