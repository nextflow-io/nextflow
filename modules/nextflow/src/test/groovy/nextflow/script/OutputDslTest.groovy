package nextflow.script

import java.nio.file.Files
import java.util.concurrent.TimeoutException

import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import spock.lang.Specification
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class OutputDslTest extends Specification {

    def 'should publish workflow outputs'() {
        given:
        def root = Files.createTempDirectory('test')
        def workDir = root.resolve('work')
        def work1 = workDir.resolve('ab/1234'); Files.createDirectories(work1)
        def work2 = workDir.resolve('cd/5678'); Files.createDirectories(work2)
        def file1 = work1.resolve('file1.txt'); file1.text = 'Hello'
        def file2 = work2.resolve('file2.txt'); file2.text = 'world'
        def target = root.resolve('results')
        and:
        def session = Mock(Session) {
            getConfig() >> [:]
            getWorkDir() >> workDir
        }
        Global.session = session
        and:
        def ch1 = new DataflowQueue()
        ch1.bind(file1)
        ch1.bind(Channel.STOP)
        and:
        def ch2 = new DataflowQueue()
        ch2.bind(file2)
        ch2.bind(Channel.STOP)
        and:
        def targets = [
            (ch1): 'foo',
            (ch2): 'bar'
        ]
        def dsl = new OutputDsl()
        and:
        SysEnv.push(NXF_FILE_ROOT: root.toString())

        when:
        dsl.directory('results')
        dsl.mode('symlink')
        dsl.overwrite(true)
        dsl.target('bar') {
            path('barbar')
            index {
                path 'index.csv'
            }
        }
        dsl.build(targets)

        def now = System.currentTimeMillis()
        while( !dsl.complete ) {
            sleep 100
            if( System.currentTimeMillis() - now > 5_000 )
                throw new TimeoutException()
        }

        then:
        target.resolve('foo/file1.txt').text == 'Hello'
        target.resolve('barbar/file2.txt').text == 'world'
        target.resolve('barbar/index.csv').text == """\
            "file2","${target}/barbar/file2.txt"
            """.stripIndent()
        and:
        1 * session.notifyFilePublish(target.resolve('foo/file1.txt'), file1)
        1 * session.notifyFilePublish(target.resolve('barbar/file2.txt'), file2)
        1 * session.notifyFilePublish(target.resolve('barbar/index.csv'))

        cleanup:
        SysEnv.pop()
        root?.deleteDir()
    }

}
