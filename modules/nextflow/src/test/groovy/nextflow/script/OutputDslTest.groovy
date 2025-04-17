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
        def outputDir = root.resolve('results')
        def workDir = root.resolve('work')
        def work1 = workDir.resolve('ab/1234'); Files.createDirectories(work1)
        def work2 = workDir.resolve('cd/5678'); Files.createDirectories(work2)
        def file1 = work1.resolve('file1.txt'); file1.text = 'Hello'
        def file2 = work2.resolve('file2.txt'); file2.text = 'world'
        and:
        def session = Mock(Session) {
            getOutputs() >> [:]
            getConfig() >> [
                workflow: [
                    output: [
                        mode: 'symlink',
                        overwrite: true
                    ]
                ]
            ]
            getOutputDir() >> outputDir
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
        session.outputs.put('foo', ch1)
        session.outputs.put('bar', ch2)
        def dsl = new OutputDsl()
        and:
        SysEnv.push(NXF_FILE_ROOT: root.toString())

        when:
        dsl.declare('foo') {
            path('foo')
        }
        dsl.declare('bar') {
            path { v -> "${'barbar'}" }
            index {
                path 'index.csv'
            }
        }
        dsl.apply(session)

        def now = System.currentTimeMillis()
        while( !dsl.complete ) {
            sleep 100
            if( System.currentTimeMillis() - now > 5_000 )
                throw new TimeoutException()
        }

        then:
        outputDir.resolve('foo/file1.txt').text == 'Hello'
        outputDir.resolve('barbar/file2.txt').text == 'world'
        outputDir.resolve('index.csv').text == """\
            "${outputDir}/barbar/file2.txt"
            """.stripIndent()
        and:
        1 * session.notifyFilePublish(outputDir.resolve('foo/file1.txt'), file1, null)
        1 * session.notifyFilePublish(outputDir.resolve('barbar/file2.txt'), file2, null)
        1 * session.notifyFilePublish(outputDir.resolve('index.csv'), null, null)

        cleanup:
        SysEnv.pop()
        root?.deleteDir()
    }

    def 'should set publish options in output declaration' () {
        when:
        def dsl1 = new OutputDsl.DeclareDsl()
        then:
        dsl1.getOptions() == [:]

        when:
        def dsl2 = new OutputDsl.DeclareDsl()
        and:
        dsl2.contentType('simple/text')
        dsl2.enabled(true)
        dsl2.ignoreErrors(true)
        dsl2.mode('someMode')
        dsl2.overwrite(true)
        dsl2.storageClass('someClass')
        dsl2.tags([foo:'1',bar:'2'])
        then:
        dsl2.getOptions() == [
            contentType:'simple/text',
            enabled: true,
            ignoreErrors: true,
            mode: 'someMode',
            overwrite: true,
            storageClass: 'someClass',
            tags: [foo:'1',bar:'2']
        ]
    }

    def 'should set index directives' () {
        when:
        def dsl1 = new OutputDsl.IndexDsl()
        then:
        dsl1.getOptions() == [:]

        when:
        def dsl2 = new OutputDsl.IndexDsl()
        and:
        dsl2.header(true)
        dsl2.path('path')
        dsl2.sep(',')
        then:
        dsl2.getOptions() == [
            header: true,
            path: 'path',
            sep: ','
        ]
    }

}
