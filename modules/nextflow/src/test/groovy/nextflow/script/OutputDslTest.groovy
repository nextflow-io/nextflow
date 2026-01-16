package nextflow.script

import java.nio.file.Files

import nextflow.Channel
import nextflow.Session
import nextflow.SysEnv
import nextflow.trace.event.FilePublishEvent
import nextflow.trace.event.WorkflowOutputEvent
import nextflow.trace.event.WorkflowPublishEvent
import spock.lang.Specification

import static test.ScriptHelper.runDataflow
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
        def config = [
            outputDir: outputDir,
            workDir: workDir,
            workflow: [
                output: [
                    mode: 'symlink',
                    overwrite: true
                ]
            ]
        ]
        and:
        SysEnv.push(NXF_FILE_ROOT: root.toString())

        when:
        def session = Spy(new Session(config))

        session.outputs.put('foo', Channel.of(file1))
        session.outputs.put('bar', Channel.of(file2))

        def dsl = new OutputDsl()
        dsl.declare('foo') {
            path('foo')
        }
        dsl.declare('bar') {
            path { v -> "${'barbar'}" }
            label 'foo'
            label 'bar'
            index {
                path 'index.csv'
            }
        }
        dsl.apply(session)
        session.fireDataflowNetwork()
        dsl.getOutput()

        then:
        outputDir.resolve('foo/file1.txt').text == 'Hello'
        outputDir.resolve('barbar/file2.txt').text == 'world'
        outputDir.resolve('index.csv').text == """\
            "${outputDir}/barbar/file2.txt"
            """.stripIndent()
        and:
        session.notifyFilePublish(new FilePublishEvent(file1, outputDir.resolve('foo/file1.txt'), null))
        session.notifyFilePublish(new FilePublishEvent(file2, outputDir.resolve('barbar/file2.txt'), ['foo', 'bar']))
        session.notifyWorkflowPublish(new WorkflowPublishEvent('foo', outputDir.resolve('foo/file1.txt')))
        session.notifyWorkflowPublish(new WorkflowPublishEvent('bar', outputDir.resolve('barbar/file2.txt')))
        session.notifyWorkflowOutput(new WorkflowOutputEvent('foo', [outputDir.resolve('foo/file1.txt')], null))
        session.notifyWorkflowOutput(new WorkflowOutputEvent('bar', [outputDir.resolve('barbar/file2.txt')], outputDir.resolve('index.csv')))
        session.notifyFilePublish(new FilePublishEvent(null, outputDir.resolve('index.csv'), ['foo', 'bar']))

        cleanup:
        SysEnv.pop()
        root?.deleteDir()
    }

    def 'should accept empty output declaration'() {
        given:
        def root = Files.createTempDirectory('test')
        def outputDir = root.resolve('results')
        def workDir = root.resolve('work')
        def work1 = workDir.resolve('ab/1234'); Files.createDirectories(work1)
        def file1 = work1.resolve('file1.txt'); file1.text = 'Hello'
        and:
        def config = [
            outputDir: outputDir,
            workDir: workDir
        ]
        and:
        SysEnv.push(NXF_FILE_ROOT: root.toString())

        when:
        def session = Spy(new Session(config))

        session.outputs.put('foo', Channel.of(file1))

        def dsl = new OutputDsl()
        dsl.declare('foo') {
        }
        dsl.apply(session)
        session.fireDataflowNetwork()
        dsl.getOutput()

        then:
        outputDir.resolve('file1.txt').text == 'Hello'
        and:
        session.notifyFilePublish(new FilePublishEvent(file1, outputDir.resolve('file1.txt'), null))
        session.notifyWorkflowPublish(new WorkflowPublishEvent('foo', outputDir.resolve('file1.txt')))
        session.notifyWorkflowOutput(new WorkflowOutputEvent('foo', [outputDir.resolve('file1.txt')], null))

        cleanup:
        SysEnv.pop()
        root?.deleteDir()
    }

    def 'should preserve non-task output files in workflow output'() {
        given:
        def root = Files.createTempDirectory('test')
        def outputDir = root.resolve('results')
        def workDir = root.resolve('work')
        def inputDir = root.resolve('inputs'); Files.createDirectories(inputDir)
        def file1 = inputDir.resolve('file1.txt'); file1.text = 'Hello'
        def file2 = inputDir.resolve('file2.txt'); file2.text = 'world'
        def record = [id: '1', file1: file1, file2: file2]
        and:
        def config = [
            outputDir: outputDir,
            workDir: workDir
        ]
        and:
        SysEnv.push(NXF_FILE_ROOT: root.toString())

        when:
        def session = Spy(new Session(config))

        session.outputs.put('foo', Channel.of(record))

        def dsl = new OutputDsl()
        dsl.declare('foo') {
        }
        dsl.apply(session)
        session.fireDataflowNetwork()
        dsl.getOutput()

        then:
        0 * session.notifyFilePublish(_)
        session.notifyWorkflowPublish(new WorkflowPublishEvent('foo', record))
        session.notifyWorkflowOutput(new WorkflowOutputEvent('foo', [ record ], null))

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
            tags: [foo:'1',bar:'2'],
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
            sep: ',',
        ]
    }

}
