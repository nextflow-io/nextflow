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
}
