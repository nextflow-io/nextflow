package nextflow

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SessionTest extends Specification {


    def testBaseDirAndBinDir() {

        setup:
        def base = File.createTempDir()
        def bin = new File(base,'bin'); bin.mkdir()

        when:
        def session = new Session()
        then:
        session.baseDir == null
        session.binDir == null

        when:
        session = new Session()
        session.baseDir = new File('some/folder')
        then:
        session.baseDir == new File('some/folder')
        session.binDir == null

        when:
        session.baseDir = base
        then:
        session.baseDir == base
        session.binDir.exists()

        cleanup:
        base.deleteDir()

    }


    def testGetQueueSize() {

        setup:
        def session = [:] as Session

        when:
        session.config = [ executor:[sge:[queueSize: 123] ] ]
        then:
        session.getQueueSize('sge', 1) == 123
        session.getQueueSize('xxx', 1) == 1
        session.getQueueSize(null, 1) == 1

        when:
        session.config = [ executor:[ queueSize: 321, sge:[queueSize:789] ] ]
        then:
        session.getQueueSize('sge', 2) == 789
        session.getQueueSize('xxx', 2) == 321
        session.getQueueSize(null, 2) == 321


        when:
        session.config = [ executor: 'sge' ]
        then:
        session.getQueueSize('sge', 1) == 1
        session.getQueueSize('xxx', 2) == 2
        session.getQueueSize(null, 3) == 3


    }

    def testGetPollInterval() {

        setup:
        def session = [:] as Session

        when:
        session.config = [ executor:[sge:[pollInterval: 345] ] ]
        then:
        session.getPollIntervalMillis('sge') == 345
        session.getPollIntervalMillis('xxx') == 1_000
        session.getPollIntervalMillis(null) == 1_000
        session.getPollIntervalMillis(null, 2_000) == 2_000

        when:
        session.config = [ executor:[ pollInterval: 321, sge:[pollInterval:789] ] ]
        then:
        session.getPollIntervalMillis('sge') == 789
        session.getPollIntervalMillis('xxx') == 321
        session.getPollIntervalMillis(null) == 321

        when:
        session.config = [ executor: 'lsf' ]
        then:
        session.getPollIntervalMillis('sge', 33 ) == 33
        session.getPollIntervalMillis('xxx', 44 ) == 44
        session.getPollIntervalMillis(null, 55 ) == 55

    }
}
