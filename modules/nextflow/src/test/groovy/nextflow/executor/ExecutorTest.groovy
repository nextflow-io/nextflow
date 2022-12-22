package nextflow.executor

import java.nio.file.Paths

import nextflow.Session
import spock.lang.Specification
import test.TestHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ExecutorTest extends Specification {


    def 'should return stage dir' () {
        given:
        def uid = UUID.randomUUID()
        def session = Mock(Session) { getUniqueId()>>uid }
        and:
        def WORK_DIR = Paths.get('/the/work/dir')
        def executor = Spy(Executor)

        when:
        executor.getWorkDir() >> WORK_DIR
        executor.getSession() >> session
        then:
        executor.getStageDir() == WORK_DIR.resolve("stage-$uid")
    }

    def 'should check foreign file' () {
        given:
        def uid = UUID.randomUUID()
        def session = Mock(Session) { getUniqueId()>>uid }
        and:
        def local = Paths.get('/work/dir')
        def foreign1 = TestHelper.createInMemTempFile('hola.txt', 'hola mundo!')
        def foreign2 = TestHelper.createInMemTempFile('ciao.txt', 'ciao mondo!')
        and:
        def executor = Spy(Executor)

        when:
        executor.getWorkDir() >> local
        executor.getSession() >> session
        then:
        !executor.isForeignFile(local.resolve('foo'))
        executor.isForeignFile(foreign1)
        executor.isForeignFile(foreign2)

    }

}
