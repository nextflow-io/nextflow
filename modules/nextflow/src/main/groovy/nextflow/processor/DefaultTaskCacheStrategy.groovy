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

import java.nio.file.Path
import java.util.concurrent.locks.Lock

import com.google.common.hash.HashCode
import com.google.common.util.concurrent.Striped
import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.file.FileHelper
import nextflow.util.HashBuilder
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * The default {@link TaskCacheStrategy}.
 *
 * <p>Every attempt of a task is keyed by the same hash -- the task hash folded with {@code 1}, as it
 * always was for a first attempt -- so a task resumes with ONE lookup: the entry under that key names
 * the work directory of the attempt that completed. The attempt is a suffix of the work directory
 * instead: {@code ab/cdef…} for the first, {@code ab/cdef…-N} after it. A failed attempt is recorded
 * under a key of its own by {@link nextflow.cache.CacheDB#putTaskAsync}, so the task key only ever
 * holds a completed attempt.
 *
 * <p>Caches written before this layout keyed attempt {@code N} by {@code H(H(h,1),…,N)} instead: a
 * failed entry under the task key is followed along that chain, so their retries resume too.
 *
 * <p>This strategy is not a plugin extension: it is the fallback {@link TaskProcessor} uses when no
 * registered {@link TaskCacheStrategy} applies to the session, so {@link #isEnabled} is always true.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class DefaultTaskCacheStrategy implements TaskCacheStrategy {

    /**
     * Deliberately {@link TaskProcessor}'s logger rather than this class's. The two lines below --
     * the {@code Cacheable folder=…} trace and the {@code Unable to resume cached task} warning --
     * were emitted from {@code TaskProcessor.checkCachedOrLaunchTask} before this strategy was split
     * out of it, and the usual support instruction for a resume problem is
     * {@code -trace nextflow.processor.TaskProcessor}. Moving them to a new logger would silently
     * drop them from every log collected that way.
     */
    private static final Logger log = LoggerFactory.getLogger(TaskProcessor)

    /**
     * Serializes the work directory check-and-create of tasks with the same hash. A fixed set
     * of stripes keeps memory bounded: distinct hashes sharing a stripe merely contend briefly
     */
    private static final Striped<Lock> locks = Striped.lock(1024)

    @Override
    boolean isEnabled(Session session) { return true }

    @Override
    void resolve(TaskRun task, HashCode hash, boolean shouldTryCache, TaskResolver resolver) {
        // a retry receives the key it was launched under, a first attempt the task hash
        final key = task.runType == TaskProcessor.RunType.RETRY
            ? hash
            : HashBuilder.defaultHasher().putBytes(hash.asBytes()).putInt(1).hash()

        if( shouldTryCache && resume(task, key, resolver) )
            return

        final workDir = claimWorkDir(resolver.workDirFor(key), task.failCount + task.rerunCount + 1)
        resolver.launch(task, key, workDir)
    }

    /**
     * Resume a task from the entry under its key, or, when that entry is a failed attempt written by
     * a previous version, from the retry chain of that version.
     */
    private static boolean resume(TaskRun task, HashCode key, TaskResolver resolver) {
        try {
            int tries = 1
            TaskEntry entry
            while( (entry = resolver.entry(key)) != null ) {
                log.trace "[${task.lazyName()}] Cacheable folder=${entry.trace.getWorkDir()} -- try=$tries; entry=$entry"
                if( entry.trace.isCompleted() ) {
                    final resumeDir = FileHelper.asPath(entry.trace.getWorkDir())
                    return resumeDir.exists() && resolver.resume(task, key, resumeDir, entry)
                }
                key = HashBuilder.defaultHasher().putBytes(key.asBytes()).putInt(++tries).hash()
            }
        }
        catch (Throwable t) {
            log.warn1("[${task.lazyName()}] Unable to resume cached task -- See log file for details", causedBy: t)
        }
        return false
    }

    /**
     * Create the first free work directory of a task, starting at {@code attempt}: an existing one
     * belongs to a previous attempt, or to an identical task instance of this run.
     */
    private static Path claimWorkDir(Path first, int attempt) {
        final lock = locks.get(first)
        lock.lock()
        try {
            while( true ) {
                final workDir = attempt == 1 ? first : first.resolveSibling("${first.fileName}-${attempt}".toString())
                if( !workDir.exists() ) {
                    if( !workDir.mkdirs() )
                        throw new IOException("Unable to create directory=$workDir -- check file system permissions")
                    return workDir
                }
                attempt++
            }
        }
        finally {
            lock.unlock()
        }
    }

}
