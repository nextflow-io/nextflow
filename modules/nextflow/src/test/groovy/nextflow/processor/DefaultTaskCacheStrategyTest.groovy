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
 * The default resolution loop, as it was in {@code TaskProcessor.checkCachedOrLaunchTask}: the
 * {@code tries} fold, the entry lookup, the resume, the bump past an existing work dir, the
 * creation of the first free one.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DefaultTaskCacheStrategyTest extends Specification {

    Path root

    def setup() { root = Files.createTempDirectory('wd') }

    def cleanup() { root?.deleteDir() }

    /**
     * The per-attempt hash the loop derives for attempt {@code tries}: the task hash folded with the
     * attempt number, CHAINED -- each attempt re-hashes the previous attempt's hash, so attempt 2 is
     * {@code H(H(hash,1),2)}, and a task retried after {@code failCount} failures starts the chain at
     * {@code H(hash, failCount+1)}.
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
        }
    }

    private TaskEntry entryIn(Path dir, boolean completed = true) {
        new TaskEntry(Mock(TraceRecord) { getWorkDir() >> dir.toUriString(); isCompleted() >> completed }, null)
    }

    def 'is always enabled: it is the fallback, not an extension'() {
        expect:
        new DefaultTaskCacheStrategy().isEnabled(Mock(Session))
        new DefaultTaskCacheStrategy().isEnabled(null)
    }

    def 'creates the work dir of the first attempt and launches there when nothing is cached'() {
        given:
        def resolver = Mock(TaskResolver)
        def task = task()
        def hash = HashCode.fromInt(100)
        def workDir = root.resolve('ab').resolve('cdef')

        when:
        new DefaultTaskCacheStrategy().resolve(task, hash, true, resolver)

        then: 'attempt 1 is looked up under the folded hash, not the task hash'
        1 * resolver.entry(attempt(hash, 1)) >> null
        1 * resolver.workDirFor(attempt(hash, 1)) >> workDir
        and: 'the directory is created here, before the launch, and nothing is resumed'
        1 * resolver.launch(task, attempt(hash, 1), workDir) >> { assert Files.isDirectory(workDir) }
        0 * resolver.resume(*_)
        and: 'the later attempts are looked up, in case attempt 1 was removed'
        DefaultTaskCacheStrategy.MAX_LOOKAHEAD * resolver.entry({ it != attempt(hash, 1) }) >> null
    }

    def 'resumes from a completed entry whose work dir exists'() {
        given:
        def resolver = Mock(TaskResolver)
        def task = task()
        def hash = HashCode.fromInt(100)
        def resumeDir = Files.createDirectories(root.resolve('ab').resolve('cdef'))
        def entry = entryIn(resumeDir)

        when:
        new DefaultTaskCacheStrategy().resolve(task, hash, true, resolver)

        then: 'the task is handed to the resolver, which resumes a copy of it (see TaskResolver.resume)'
        1 * resolver.entry(attempt(hash, 1)) >> entry
        1 * resolver.resume(task, attempt(hash, 1), resumeDir, entry) >> true
        0 * resolver.launch(*_)
        0 * resolver.workDirFor(_)
    }

    def 'bumps past an attempt whose work dir exists: not resuming, a stale entry, or a failed resume'() {
        given:
        def resolver = Mock(TaskResolver)
        def task = task()
        def hash = HashCode.fromInt(100)
        def usedDir = Files.createDirectories(root.resolve('ab').resolve('cdef'))
        def freeDir = root.resolve('gh').resolve('ijkl')

        when: 'the run is not resuming'
        new DefaultTaskCacheStrategy().resolve(task, hash, false, resolver)
        then: 'the existing entry is skipped without a resume and the next attempt is launched'
        1 * resolver.entry(attempt(hash, 1)) >> entryIn(usedDir)
        0 * resolver.resume(*_)
        1 * resolver.entry(attempt(hash, 2)) >> null
        1 * resolver.workDirFor(attempt(hash, 2)) >> freeDir
        1 * resolver.launch(task, attempt(hash, 2), freeDir)

        when: 'the entry did not complete'
        freeDir.deleteDir()
        new DefaultTaskCacheStrategy().resolve(task, hash, true, resolver)
        then:
        1 * resolver.entry(attempt(hash, 1)) >> entryIn(usedDir, false)
        0 * resolver.resume(*_)
        1 * resolver.entry(attempt(hash, 2)) >> null
        1 * resolver.workDirFor(attempt(hash, 2)) >> freeDir
        1 * resolver.launch(task, attempt(hash, 2), freeDir)

        when: 'the entry completed but its outputs cannot be resumed'
        freeDir.deleteDir()
        new DefaultTaskCacheStrategy().resolve(task, hash, true, resolver)
        then:
        1 * resolver.entry(attempt(hash, 1)) >> entryIn(usedDir)
        1 * resolver.resume(task, attempt(hash, 1), usedDir, _) >> false
        1 * resolver.entry(attempt(hash, 2)) >> null
        1 * resolver.workDirFor(attempt(hash, 2)) >> freeDir
        1 * resolver.launch(task, attempt(hash, 2), freeDir)

        when: 'the resume blows up'
        freeDir.deleteDir()
        new DefaultTaskCacheStrategy().resolve(task, hash, true, resolver)
        then: 'it is a warning, and the task re-executes in the next attempt'
        1 * resolver.entry(attempt(hash, 1)) >> entryIn(usedDir)
        1 * resolver.resume(task, attempt(hash, 1), usedDir, _) >> { throw new IOException('transient 500') }
        1 * resolver.entry(attempt(hash, 2)) >> null
        1 * resolver.workDirFor(attempt(hash, 2)) >> freeDir
        1 * resolver.launch(task, attempt(hash, 2), freeDir)
        noExceptionThrown()
    }

    def 'bumps past a work dir that exists without an entry -- an identical instance of this run'() {
        given:
        def resolver = Mock(TaskResolver)
        def task = task()
        def hash = HashCode.fromInt(100)
        def usedDir = Files.createDirectories(root.resolve('ab').resolve('cdef'))
        def freeDir = root.resolve('gh').resolve('ijkl')

        when:
        new DefaultTaskCacheStrategy().resolve(task, hash, true, resolver)

        then:
        1 * resolver.entry(attempt(hash, 1)) >> null
        1 * resolver.workDirFor(attempt(hash, 1)) >> usedDir
        1 * resolver.workDirFor(attempt(hash, 2)) >> freeDir
        1 * resolver.launch(task, attempt(hash, 2), freeDir)
        and: 'looked up by the walk and by the look-ahead'
        _ * resolver.entry(_) >> null
    }

    def 'starts counting attempts after the failures of a retried task'() {
        given:
        def resolver = Mock(TaskResolver)
        def task = task(2)
        def hash = HashCode.fromInt(100)
        def workDir = root.resolve('ab').resolve('cdef')

        when:
        new DefaultTaskCacheStrategy().resolve(task, hash, false, resolver)

        then:
        1 * resolver.entry(attempt(hash, 3, 2)) >> null
        1 * resolver.workDirFor(attempt(hash, 3, 2)) >> workDir
        1 * resolver.launch(task, attempt(hash, 3, 2), workDir)
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
        1 * resolver.entry(_) >> null
        1 * resolver.workDirFor(_) >> blocker.resolve('xy')
        0 * resolver.launch(*_)
        def err = thrown(IOException)
        err.message.startsWith('Unable to create directory=')
    }

    def 'steps past a missing attempt to a later cached one -- e.g. a failed attempt removed by clean'() {
        given:
        def resolver = Mock(TaskResolver)
        def task = task()
        def hash = HashCode.fromInt(100)
        def resumeDir = Files.createDirectories(root.resolve('ab').resolve('cdef'))
        def entry = entryIn(resumeDir)

        when:
        new DefaultTaskCacheStrategy().resolve(task, hash, true, resolver)

        then: 'attempts 1 and 2 are gone, attempt 3 is cached'
        resolver.entry(attempt(hash, 3)) >> entry
        resolver.entry(_) >> null
        1 * resolver.resume(task, attempt(hash, 3), resumeDir, entry) >> true
        0 * resolver.launch(*_)
    }

    def 'does not look ahead further than the limit, nor when not resuming'() {
        given:
        def resolver = Mock(TaskResolver)
        def hash = HashCode.fromInt(100)
        def beyond = attempt(hash, DefaultTaskCacheStrategy.MAX_LOOKAHEAD + 2)
        def entry = entryIn(Files.createDirectories(root.resolve('ab').resolve('cdef')))
        def workDir = root.resolve('gh').resolve('ijkl')

        when:
        new DefaultTaskCacheStrategy().resolve(task(), hash, true, resolver)
        then:
        resolver.entry(beyond) >> entry
        resolver.entry(_) >> null
        1 * resolver.workDirFor(attempt(hash, 1)) >> workDir
        1 * resolver.launch(_, attempt(hash, 1), workDir)
        0 * resolver.resume(*_)

        when:
        workDir.deleteDir()
        new DefaultTaskCacheStrategy().resolve(task(), hash, false, resolver)
        then:
        1 * resolver.entry(attempt(hash, 1)) >> null
        0 * resolver.entry(_)
        1 * resolver.workDirFor(attempt(hash, 1)) >> workDir
        1 * resolver.launch(_, attempt(hash, 1), workDir)
    }

    def 'a retry continues the walk of a task launched past an existing work dir'() {
        given:
        def strategy = new DefaultTaskCacheStrategy()
        def resolver = Mock(TaskResolver)
        def hash = HashCode.fromInt(100)
        def leftover = Files.createDirectories(root.resolve('ab').resolve('cdef'))
        def first = root.resolve('gh').resolve('ijkl')
        def retry = root.resolve('mn').resolve('opqr')

        when: 'attempt 1 is a leftover of a previous run, so the task launches at attempt 2'
        strategy.resolve(task(), hash, true, resolver)
        then:
        resolver.entry(attempt(hash, 1)) >> entryIn(leftover, false)
        resolver.entry(_) >> null
        1 * resolver.workDirFor(attempt(hash, 2)) >> first
        1 * resolver.launch(_, attempt(hash, 2), first)

        when: 'it fails: the retry is keyed at attempt 3 of the walk, where the next resume looks for it'
        strategy.resolve(task(1), attempt(hash, 2), false, resolver)
        then:
        1 * resolver.entry(attempt(hash, 3)) >> null
        1 * resolver.workDirFor(attempt(hash, 3)) >> retry
        1 * resolver.launch(_, attempt(hash, 3), retry)
    }

    def 'resumes a retry written off the walk by a previous version'() {
        given:
        def resolver = Mock(TaskResolver)
        def hash = HashCode.fromInt(100)
        def failed1 = Files.createDirectories(root.resolve('ab').resolve('0001'))
        def failed2 = Files.createDirectories(root.resolve('ab').resolve('0002'))
        def resumeDir = Files.createDirectories(root.resolve('ab').resolve('0003'))
        def entry = entryIn(resumeDir)
        and: 'launched at attempt 2 of the walk, then retried at H(attempt 2, 2)'
        def legacy = HashBuilder.defaultHasher().putBytes(attempt(hash, 2).asBytes()).putInt(2).hash()

        when:
        new DefaultTaskCacheStrategy().resolve(task(), hash, true, resolver)

        then:
        resolver.entry(attempt(hash, 1)) >> entryIn(failed1, false)
        resolver.entry(attempt(hash, 2)) >> entryIn(failed2, false)
        resolver.entry(legacy) >> entry
        1 * resolver.resume(_, legacy, resumeDir, entry) >> true
        0 * resolver.launch(*_)
    }

}
