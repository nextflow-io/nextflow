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

import java.nio.file.CopyOption
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.TimeUnit

import nextflow.Session
import nextflow.exception.ScriptRuntimeException
import nextflow.trace.event.FilePublishEvent
import spock.lang.Specification
import test.TestHelper
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
        def publisher = new FilePublisher(session, [mode: 'copy'])

        when:
        publisher.publish(file1, outputDir.resolve('alpha.txt'))
        publisher.publish(file2, outputDir.resolve('beta.txt'))
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
        def publisher = new FilePublisher(session, [mode: 'copy'])

        when:
        publisher.publish(file1, outputDir.resolve('copied.txt'))
        awaitPublish()

        then:
        outputDir.resolve('copied.txt').text == 'Hello'

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
        def publisher = new FilePublisher(session, [mode: 'copy'])

        when:
        publisher.publish(file1, outputDir.resolve('one.txt'))
        publisher.publish(file1, outputDir.resolve('two.txt'))
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
        def publisher = new FilePublisher(session, [:])

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
        new FilePublisher(session, [mode: 'nope'])

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.startsWith "Invalid publish mode 'nope'"
    }

    def 'should resolve the publish mode'() {
        given:
        def session = mockSession()
        def source = Path.of('/work/ab/1234/file.txt')
        def target = Path.of('/results/file.txt')

        expect:
        new FilePublisher(session, OPTS).resolveMode(source, target) == EXPECTED

        where:
        OPTS                        | EXPECTED
        [:]                         | FilePublisher.Mode.SYMLINK
        [defaultMode: 'rellink']    | FilePublisher.Mode.RELLINK
        [mode: 'copy']              | FilePublisher.Mode.COPY
        [mode: 'move']              | FilePublisher.Mode.MOVE
        [mode: 'rellink']           | FilePublisher.Mode.RELLINK
        [mode: 'copyNoFollow']      | FilePublisher.Mode.COPY_NO_FOLLOW
    }

    def 'should detect a publish mode mismatch'() {
        given:
        def folder = Files.createTempDirectory('test')
        def realFile = folder.resolve('real.txt'); realFile.text = 'Hello'
        def symlink = folder.resolve('link.txt')
        Files.createSymbolicLink(symlink, realFile)
        def target = TARGET == 'symlink' ? symlink : realFile

        expect:
        new FilePublisher(mockSession(), [:]).checkPublishModeMismatch(target, FilePublisher.Mode.valueOf(MODE)) == EXPECTED

        cleanup:
        folder?.deleteDir()

        where:
        MODE        | TARGET    | EXPECTED
        'COPY'      | 'symlink' | true
        'MOVE'      | 'symlink' | true
        'SYMLINK'   | 'file'    | true
        'RELLINK'   | 'file'    | true
        'COPY'      | 'file'    | false
        'SYMLINK'   | 'symlink' | false
        'RELLINK'   | 'symlink' | false
        // a hard link is a regular file, not a symlink
        'LINK'      | 'file'    | false
        'LINK'      | 'symlink' | true
    }

    def 'should change mode to copy when the target is a foreign file system'() {
        given:
        def source = TestHelper.createInMemTempDir().resolve('file.txt')
        def target = TestHelper.createInMemTempDir().resolve('file.txt')

        expect:
        new FilePublisher(mockSession(), OPTS).resolveMode(source, target) == FilePublisher.Mode.COPY

        where:
        OPTS << [ [:], [mode: 'symlink'], [mode: 'link'], [mode: 'rellink'], [defaultMode: 'rellink'] ]
    }

    def 'should check same real path'() {
        given:
        def folder = Files.createTempDirectory('test')
        def pubDir = folder.resolve('pub-dir'); pubDir.mkdir()
        def workDir = folder.resolve('work-dir'); workDir.mkdir()
        def foo = workDir.resolve('foo.txt'); foo.text = 'This is foo'
        def bar = workDir.resolve('bar.txt'); bar.text = 'This is bar'
        def linkToBar = Files.createSymbolicLink(pubDir.resolve('link-to-bar'), bar)
        def publisher = new FilePublisher(mockSession(), [:])

        expect:
        publisher.checkIsSameRealPath(bar, linkToBar, FilePublisher.Mode.SYMLINK)
        !publisher.checkIsSameRealPath(bar, foo, FilePublisher.Mode.SYMLINK)
        !publisher.checkIsSameRealPath(bar, linkToBar, FilePublisher.Mode.COPY)

        cleanup:
        folder?.deleteDir()
    }

    def 'should detect targets that overlap with the source directory'() {
        given:
        def folder = Files.createTempDirectory('test')
        def pubDir = folder.resolve('pub-dir'); pubDir.mkdir()
        def workDir = folder.resolve('work-dir'); workDir.mkdir()
        def publisher = new FilePublisher(mockSession(), [sourceDir: workDir])
        def symlink = FilePublisher.Mode.SYMLINK

        when:
        def foo = pubDir.resolve('foo.txt'); foo.text = 'This is foo'
        def bar = workDir.resolve('bar.txt'); bar.text = 'This is bar'
        then:
        !publisher.checkSourcePathConflicts(foo, symlink)
        publisher.checkSourcePathConflicts(bar, symlink)
        !publisher.checkSourcePathConflicts(bar, FilePublisher.Mode.COPY)
        !new FilePublisher(mockSession(), [:]).checkSourcePathConflicts(bar, symlink)

        when:
        def linkOK = Files.createSymbolicLink(pubDir.resolve('link1.txt'), foo)
        then:
        !publisher.checkSourcePathConflicts(linkOK, symlink)

        when:
        def linkNotOK = Files.createSymbolicLink(pubDir.resolve('link2.txt'), bar)
        then:
        publisher.checkSourcePathConflicts(linkNotOK, symlink)

        cleanup:
        folder?.deleteDir()
    }

    def 'should re-publish when the publish mode of an existing target changes'() {
        given:
        def folder = Files.createTempDirectory('test')
        def pubDir = folder.resolve('pub-dir'); pubDir.mkdir()
        def workDir = folder.resolve('work-dir'); workDir.mkdir()
        def source = workDir.resolve('foo.txt'); source.text = 'Hello'
        def target = pubDir.resolve('foo.txt')
        def session = mockSession()
        def opts = [sourceDir: workDir, overwrite: 'standard']

        when:
        new FilePublisher(session, opts).publishFile(source, target, FilePublisher.Mode.RELLINK)
        then:
        Files.isSymbolicLink(target)

        when:
        new FilePublisher(session, opts).publishFile(source, target, FilePublisher.Mode.COPY)
        then:
        !Files.isSymbolicLink(target)
        target.text == 'Hello'

        when:
        new FilePublisher(session, opts).publishFile(source, target, FilePublisher.Mode.RELLINK)
        then:
        Files.isSymbolicLink(target)

        when: 'an explicit `overwrite false` is honored even on a mode mismatch'
        new FilePublisher(session, opts + [overwrite: false]).publishFile(source, target, FilePublisher.Mode.COPY)
        then:
        Files.isSymbolicLink(target)

        cleanup:
        folder?.deleteDir()
    }

    def 'should return copy options'() {
        given:
        def session = Mock(Session) { getConfig() >> CONFIG }
        def publisher = new FilePublisher(session, [:])

        expect:
        publisher.copyOpts() == EXPECTED as CopyOption[]
        publisher.copyOpts(LinkOption.NOFOLLOW_LINKS) == ([LinkOption.NOFOLLOW_LINKS] + EXPECTED) as CopyOption[]

        where:
        CONFIG                                          | EXPECTED
        [:]                                             | []
        [workflow: [output: [copyAttributes: true]]]    | [StandardCopyOption.COPY_ATTRIBUTES]
    }

}
