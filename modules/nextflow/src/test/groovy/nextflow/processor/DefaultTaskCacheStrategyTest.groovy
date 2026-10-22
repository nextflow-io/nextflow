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

package nextflow.processor

import java.nio.file.Files
import java.nio.file.Path

import com.google.common.hash.HashCode
import nextflow.Session
import nextflow.trace.TraceRecord
import nextflow.util.HashBuilder
import spock.lang.Specification

/**
 * The default resolution: one key for every attempt of a task, the resume from the work dir its entry
 * names, the attempt suffix of the work dir, and the retry chain of caches written by previous versions.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DefaultTaskCacheStrategyTest extends Specification {

    Path root

    def setup() { root = Files.createTempDirectory('wd') }

    def cleanup() { root?.deleteDir() }

    /**
     * The per-attempt hash of previous versions: the task hash folded with the attempt number, CHAINED,
     * so attempt 2 is {@code H(H(hash,1),2)}. Attempt 1 is the key of every attempt now.
     */
    private static HashCode attempt(HashCode hash, int tries, int failCount = 0) {
        def h = hash
        for( int t = failCount + 1; t <= tries; t++ )
            h = HashBuilder.defaultHasher().putBytes(h.asBytes()).putInt(t).hash()
        return h
    }

    private TaskRun task(int failCount = 0) {
        Mock(TaskRun) {
            getFailCount() >> failCount
            getRunType() >> (failCount ? TaskProcessor.RunType.RETRY : TaskProcessor.RunType.SUBMIT)
        }
    }

    private TaskEntry entryIn(Path dir, boolean completed = true) {
        new TaskEntry(Mock(TraceRecord) { getWorkDir() >> dir.toUriString(); isCompleted() >> completed }, null)
    }

    def 'is always enabled: it is the fallback, not an extension'() {
        expect:
        new DefaultTaskCacheStrategy().isEnabled(null)
        new DefaultTaskCacheStrategy().isEnabled(Mock(Session))
    }

    def 'creates the work dir of the first attempt and launches there when nothing is cached'() {
        given:
        def resolver = Mock(TaskResolver)
        def task = task()
        def hash = HashCode.fromInt(100)
        def workDir = root.resolve('ab').resolve('cdef')

        when:
        new DefaultTaskCacheStrategy().resolve(task, hash, true, resolver)

        then: 'the task is keyed by the hash folded with 1, as a first attempt always was'
        1 * resolver.entry(attempt(hash, 1)) >> null
        1 * resolver.workDirFor(attempt(hash, 1)) >> workDir
        and: 'the directory is created here, before the launch, and nothing is resumed'
        1 * resolver.launch(task, attempt(hash, 1), workDir) >> { assert Files.isDirectory(workDir) }
        0 * resolver.resume(*_)
        0 * resolver.entry(_)
    }

    def 'resumes from the work dir named by the entry, whichever attempt completed'() {
        given:
        def resolver = Mock(TaskResolver)
        def task = task()
        def hash = HashCode.fromInt(100)
        def resumeDir = Files.createDirectories(root.resolve('ab').resolve('cdef-3'))
        def entry = entryIn(resumeDir)

        when:
        new DefaultTaskCacheStrategy().resolve(task, hash, true, resolver)

        then: 'ONE lookup -- the failed attempts are not needed'
        1 * resolver.entry(attempt(hash, 1)) >> entry
        1 * resolver.resume(task, attempt(hash, 1), resumeDir, entry) >> true
        0 * resolver.entry(_)
        0 * resolver.launch(*_)
        0 * resolver.workDirFor(_)
    }

    def 'launches in the next free work dir when the entry cannot be resumed'() {
        given:
        def resolver = Mock(TaskResolver)
        def task = task()
        def hash = HashCode.fromInt(100)
        def first = Files.createDirectories(root.resolve('ab').resolve('cdef'))
        def next = root.resolve('ab').resolve('cdef-2')

        when: 'the run is not resuming'
        new DefaultTaskCacheStrategy().resolve(task, hash, false, resolver)
        then: 'no lookup, and the existing directory is skipped'
        0 * resolver.entry(_)
        1 * resolver.workDirFor(attempt(hash, 1)) >> first
        1 * resolver.launch(task, attempt(hash, 1), next)

        when: 'the completed entry names a work dir that is gone'
        next.deleteDir()
        new DefaultTaskCacheStrategy().resolve(task, hash, true, resolver)
        then:
        1 * resolver.entry(attempt(hash, 1)) >> entryIn(root.resolve('gone'))
        0 * resolver.resume(*_)
        1 * resolver.workDirFor(attempt(hash, 1)) >> first
        1 * resolver.launch(task, attempt(hash, 1), next)

        when: 'its outputs cannot be resumed'
        next.deleteDir()
        new DefaultTaskCacheStrategy().resolve(task, hash, true, resolver)
        then:
        1 * resolver.entry(attempt(hash, 1)) >> entryIn(first)
        1 * resolver.resume(task, attempt(hash, 1), first, _) >> false
        1 * resolver.workDirFor(attempt(hash, 1)) >> first
        1 * resolver.launch(task, attempt(hash, 1), next)

        when: 'the resume blows up'
        next.deleteDir()
        new DefaultTaskCacheStrategy().resolve(task, hash, true, resolver)
        then: 'it is a warning, and the task re-executes'
        1 * resolver.entry(attempt(hash, 1)) >> entryIn(first)
        1 * resolver.resume(task, attempt(hash, 1), first, _) >> { throw new IOException('transient 500') }
        1 * resolver.workDirFor(attempt(hash, 1)) >> first
        1 * resolver.launch(task, attempt(hash, 1), next)
        noExceptionThrown()
    }

    def 'a retry keeps the key of the task and starts at the work dir of its attempt'() {
        given:
        def resolver = Mock(TaskResolver)
        def key = attempt(HashCode.fromInt(100), 1)
        def first = root.resolve('ab').resolve('cdef')

        when: 'retried after two failures: the key passed in is the task key, not folded again'
        new DefaultTaskCacheStrategy().resolve(task(2), key, false, resolver)
        then:
        1 * resolver.workDirFor(key) >> first
        1 * resolver.launch(_, key, root.resolve('ab').resolve('cdef-3'))
        0 * resolver.entry(_)

        when: 'that work dir is in use -- e.g. by an identical task instance of this run'
        Files.createDirectories(root.resolve('ab').resolve('cdef-4'))
        new DefaultTaskCacheStrategy().resolve(task(3), key, false, resolver)
        then:
        1 * resolver.workDirFor(key) >> first
        1 * resolver.launch(_, key, root.resolve('ab').resolve('cdef-5'))
    }

    def 'resumes a retry from a cache written by a previous version'() {
        given:
        def resolver = Mock(TaskResolver)
        def task = task()
        def hash = HashCode.fromInt(100)
        def resumeDir = Files.createDirectories(root.resolve('ij').resolve('klmn'))
        def entry = entryIn(resumeDir)

        when:
        new DefaultTaskCacheStrategy().resolve(task, hash, true, resolver)

        then: 'a failed entry under the task key is followed along the old per-attempt chain'
        1 * resolver.entry(attempt(hash, 1)) >> entryIn(root.resolve('ab').resolve('cdef'), false)
        1 * resolver.entry(attempt(hash, 2)) >> entryIn(root.resolve('ef').resolve('ghij'), false)
        1 * resolver.entry(attempt(hash, 3)) >> entry
        1 * resolver.resume(task, attempt(hash, 3), resumeDir, entry) >> true
        0 * resolver.launch(*_)
    }

    def 'launches when the old per-attempt chain has no completed entry'() {
        given:
        def resolver = Mock(TaskResolver)
        def task = task()
        def hash = HashCode.fromInt(100)
        def first = Files.createDirectories(root.resolve('ab').resolve('cdef'))

        when:
        new DefaultTaskCacheStrategy().resolve(task, hash, true, resolver)

        then:
        1 * resolver.entry(attempt(hash, 1)) >> entryIn(first, false)
        1 * resolver.entry(attempt(hash, 2)) >> null
        1 * resolver.workDirFor(attempt(hash, 1)) >> first
        1 * resolver.launch(task, attempt(hash, 1), root.resolve('ab').resolve('cdef-2'))
    }

    def 'fails when the work dir cannot be created'() {
        given:
        def resolver = Mock(TaskResolver)
        def task = task()
        def hash = HashCode.fromInt(100)
        def blocker = root.resolve('blocker'); blocker.text = 'not a directory'

        when:
        new DefaultTaskCacheStrategy().resolve(task, hash, false, resolver)

        then:
        1 * resolver.workDirFor(_) >> blocker.resolve('xy')
        0 * resolver.launch(*_)
        def err = thrown(IOException)
        err.message.startsWith('Unable to create directory=')
    }

}
