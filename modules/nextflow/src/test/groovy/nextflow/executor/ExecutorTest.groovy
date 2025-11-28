package nextflow.executor

import java.nio.file.Paths

import nextflow.Session
import nextflow.SysEnv
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

    def 'should check stage file enabled' () {
        given:
        SysEnv.push(ENV)

        expect:
        Spy(EXECUTOR).isStageFileEnabled() == EXPECTED

        cleanup:
        SysEnv.pop()

        where:
        EXECUTOR             | ENV                                       | EXPECTED
        Executor             | [:]                                       | false      // base Executor defaults to false
        AbstractGridExecutor | [:]                                       | true       // grid executor defaults to true
        Executor             | [NXF_WRAPPER_STAGE_FILE_ENABLED: 'true']  | true       // env var overrides
        AbstractGridExecutor | [NXF_WRAPPER_STAGE_FILE_ENABLED: 'false'] | false      // env var overrides grid default
    }

}
