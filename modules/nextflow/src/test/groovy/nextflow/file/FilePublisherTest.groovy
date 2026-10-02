/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.file

import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.TimeUnit

import nextflow.Session
import nextflow.exception.ScriptRuntimeException
import nextflow.trace.event.FilePublishEvent
import spock.lang.Specification
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class FilePublisherTest extends Specification {

    ExecutorService pool

    def cleanup() {
        pool?.shutdownNow()
    }

    Session mockSession() {
        pool = Executors.newSingleThreadExecutor()
        return Mock(Session) {
            getConfig() >> [:]
            publishDirExecutorService() >> pool
        }
    }

    /**
     * Wait for the publishes handed to the thread pool to complete.
     */
    void awaitPublish() {
        pool.shutdown()
        pool.awaitTermination(30, TimeUnit.SECONDS)
    }

    def 'should publish two files with the same name from different directories'() {
        given:
        def root = Files.createTempDirectory('test')
        def work1 = root.resolve('work/ab/1234'); Files.createDirectories(work1)
        def work2 = root.resolve('work/cd/5678'); Files.createDirectories(work2)
        def file1 = work1.resolve('report.txt'); file1.text = 'Hello'
        def file2 = work2.resolve('report.txt'); file2.text = 'world'
        def outputDir = root.resolve('results')
        and:
        def session = mockSession()
        def publisher = new FilePublisher(session, 'foo', [mode: 'copy'])

        when:
        publisher.publish([
            (file1): outputDir.resolve('alpha.txt'),
            (file2): outputDir.resolve('beta.txt')
        ])
        awaitPublish()

        then:
        outputDir.resolve('alpha.txt').text == 'Hello'
        outputDir.resolve('beta.txt').text == 'world'
        and:
        1 * session.notifyFilePublish(new FilePublishEvent(file1, outputDir.resolve('alpha.txt'), null))
        1 * session.notifyFilePublish(new FilePublishEvent(file2, outputDir.resolve('beta.txt'), null))

        cleanup:
        root?.deleteDir()
    }

    def 'should publish a file that does not reside in the work directory'() {
        given:
        def root = Files.createTempDirectory('test')
        def inputDir = root.resolve('inputs'); Files.createDirectories(inputDir)
        def file1 = inputDir.resolve('input.txt'); file1.text = 'Hello'
        def outputDir = root.resolve('results')
        and:
        def session = mockSession()
        def publisher = new FilePublisher(session, 'foo', [mode: 'copy'])

        when:
        publisher.publish(file1, outputDir.resolve('copied.txt'))
        awaitPublish()

        then:
        outputDir.resolve('copied.txt').text == 'Hello'

        cleanup:
        root?.deleteDir()
    }

    def 'should report two files published to the same target'() {
        given:
        def root = Files.createTempDirectory('test')
        def work1 = root.resolve('work/ab/1234'); Files.createDirectories(work1)
        def work2 = root.resolve('work/cd/5678'); Files.createDirectories(work2)
        def file1 = work1.resolve('report.txt'); file1.text = 'Hello'
        def file2 = work2.resolve('report.txt'); file2.text = 'world'
        def outputDir = root.resolve('results')
        def target = outputDir.resolve('report.txt')
        and:
        def session = mockSession()
        def publisher = new FilePublisher(session, 'foo', [mode: 'copy'])

        when:
        publisher.publish([(file1): target, (file2): target])
        awaitPublish()

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains "Publish target '${target.toUriString()}' for workflow output 'foo' is used by more than one file"
        and: 'nothing is published when a conflict is detected'
        !Files.exists(target)

        cleanup:
        root?.deleteDir()
    }

    def 'should publish the same file to more than one target'() {
        given:
        def root = Files.createTempDirectory('test')
        def work1 = root.resolve('work/ab/1234'); Files.createDirectories(work1)
        def file1 = work1.resolve('report.txt'); file1.text = 'Hello'
        def outputDir = root.resolve('results')
        and:
        def session = mockSession()
        def publisher = new FilePublisher(session, 'foo', [mode: 'copy'])

        when:
        publisher.publish([(file1): outputDir.resolve('one.txt')])
        publisher.publish([(file1): outputDir.resolve('two.txt')])
        awaitPublish()

        then:
        outputDir.resolve('one.txt').text == 'Hello'
        outputDir.resolve('two.txt').text == 'Hello'

        cleanup:
        root?.deleteDir()
    }

    def 'should create a symlink by default'() {
        given:
        def root = Files.createTempDirectory('test')
        def work1 = root.resolve('work/ab/1234'); Files.createDirectories(work1)
        def file1 = work1.resolve('report.txt'); file1.text = 'Hello'
        def outputDir = root.resolve('results')
        and:
        def session = mockSession()
        def publisher = new FilePublisher(session, 'foo', [:])

        when:
        publisher.publish(file1, outputDir.resolve('report.txt'))
        awaitPublish()

        then:
        Files.isSymbolicLink(outputDir.resolve('report.txt'))
        outputDir.resolve('report.txt').text == 'Hello'

        cleanup:
        root?.deleteDir()
    }

    def 'should reject an invalid publish mode'() {
        given:
        def session = mockSession()

        when:
        new FilePublisher(session, 'foo', [mode: 'nope'])

        then:
        def e = thrown(ScriptRuntimeException)
        e.message == "Invalid publish mode 'nope' for workflow output 'foo'"
    }

    def 'should resolve the publish mode'() {
        given:
        def session = mockSession()
        def source = Path.of('/work/ab/1234/file.txt')
        def target = Path.of('/results/file.txt')

        expect:
        new FilePublisher(session, 'foo', OPTS).resolveMode(source, target) == EXPECTED

        where:
        OPTS                    | EXPECTED
        [:]                     | FilePublisher.Mode.SYMLINK
        [mode: 'copy']          | FilePublisher.Mode.COPY
        [mode: 'move']          | FilePublisher.Mode.MOVE
        [mode: 'rellink']       | FilePublisher.Mode.RELLINK
        [mode: 'copyNoFollow']  | FilePublisher.Mode.COPY_NO_FOLLOW
    }

}
