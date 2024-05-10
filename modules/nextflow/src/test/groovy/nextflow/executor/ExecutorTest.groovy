package nextflow.executor

import java.nio.file.Path

import nextflow.Session
import nextflow.file.http.XPath
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import spock.lang.Specification
import spock.lang.Unroll
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
        def WORK_DIR = Path.of('/the/work/dir')
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
        def local = Path.of('/work/dir')
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

    @Unroll
    def 'should get array work dir' () {
        given:
        def task = Mock(TaskRun) { getWorkDir() >> WORK_DIR }
        def handler = Mock(TaskHandler) { getTask() >> task }
        def executor = Spy(Executor) { isFusionEnabled() >> FUSION }

        expect:
        executor.getChildWorkDir(handler) == EXPECTED
        
        where:
        FUSION  | WORK_DIR                           |  EXPECTED
        false   | Path.of('/work/dir')               | '/work/dir'
        false   | XPath.get('http://foo.com/work')   | 'http://foo.com/work'
        true    | XPath.get('http://foo.com/work')   | '/fusion/http/foo.com/work'
    }

    @Unroll
    def 'should get array launch command' () {
        given:
        def executor = Spy(Executor) {
            isFusionEnabled() >> FUSION
            isWorkDirDefaultFS() >> !FUSION
        }

        expect:
        executor.getChildLaunchCommand(WORK_DIR) == EXPECTED

        where:
        FUSION      | WORK_DIR              | EXPECTED
        false       | '/work/dir'           | 'bash /work/dir/.command.run 2>&1 > /work/dir/.command.log'
        true        | '/fusion/work/dir'    | 'bash /fusion/work/dir/.command.run'
    }

}
