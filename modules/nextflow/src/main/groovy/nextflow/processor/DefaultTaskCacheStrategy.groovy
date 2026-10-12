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
import java.util.concurrent.ConcurrentHashMap
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
 * The default {@link TaskCacheStrategy}: the per-run, local resolution loop as it has always been in
 * {@code TaskProcessor.checkCachedOrLaunchTask}, relocated here unchanged.
 *
 * <p>The task hash is folded with the attempt number ({@code tries}) into the per-attempt hash that
 * names the work directory and keys the cache entry. For each attempt the entry is looked up and, when
 * resuming, the task is resumed from it; an attempt whose work directory already exists -- a previous
 * failure, or an identical task instance of this run -- bumps to the next one; the first free work
 * directory is created under an in-process lock and the task is launched there.
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

    /**
     * How far the resume walk looks past a missing attempt, e.g. one removed by {@code nextflow clean},
     * for a later attempt that is still cached.
     */
    static final int MAX_LOOKAHEAD = 10

    /**
     * The walk position of tasks launched past an existing work dir, keyed by the hash they were launched
     * under, so that their retry continues the walk instead of restarting at {@code failCount + 1}.
     * Tasks launched at {@code failCount + 1} are not recorded.
     */
    private final Map<HashCode,Integer> positions = new ConcurrentHashMap<>()

    @Override
    boolean isEnabled(Session session) { return true }

    @Override
    void resolve(TaskRun task, HashCode hash, boolean shouldTryCache, TaskResolver resolver) {

        final position = positions.remove(hash)
        int tries = position != null ? position + 1 : task.failCount + 1
        while( true ) {
            hash = HashBuilder.defaultHasher().putBytes(hash.asBytes()).putInt(tries).hash()

            Path resumeDir = null
            boolean exists = false
            try {
                final entry = resolver.entry(hash)
                resumeDir = entry ? FileHelper.asPath(entry.trace.getWorkDir()) : null
                if( resumeDir )
                    exists = resumeDir.exists()

                log.trace "[${task.lazyName()}] Cacheable folder=${resumeDir?.toUriString()} -- exists=$exists; try=$tries; shouldTryCache=$shouldTryCache; entry=$entry"
                final cached = shouldTryCache && exists && entry.trace.isCompleted() && resolver.resume(task, hash, resumeDir, entry)
                if( cached )
                    break
            }
            catch (Throwable t) {
                log.warn1("[${task.lazyName()}] Unable to resume cached task -- See log file for details", causedBy: t)
            }

            if( exists ) {
                if( shouldTryCache && tries > 1 && resumeLegacyRetry(task, hash, resolver) )
                    break
                tries++
                continue
            }

            if( shouldTryCache && hasLaterAttempt(hash, tries, resolver) ) {
                tries++
                continue
            }

            final workDir = resolver.workDirFor(hash)
            final lock = locks.get(hash)
            lock.lock()
            try {
                if( resumeDir != workDir )
                    exists = workDir.exists()
                if( exists ) {
                    tries++
                    continue
                }
                else if( !workDir.mkdirs() )
                    throw new IOException("Unable to create directory=$workDir -- check file system permissions")
            }
            finally {
                lock.unlock()
            }

            if( tries != task.failCount + 1 )
                positions.put(hash, tries)

            // submit task for execution
            resolver.launch( task, hash, workDir )
            break
        }

    }

    /**
     * Whether a later attempt of the walk is still cached, when the current one is missing -- e.g. a
     * failed attempt removed by {@code nextflow clean} while the successful one was kept.
     */
    private static boolean hasLaterAttempt(HashCode hash, int tries, TaskResolver resolver) {
        for( int i = tries + 1; i <= tries + MAX_LOOKAHEAD; i++ ) {
            hash = HashBuilder.defaultHasher().putBytes(hash.asBytes()).putInt(i).hash()
            final entry = resolver.entry(hash)
            if( entry && FileHelper.asPath(entry.trace.getWorkDir()).exists() )
                return true
        }
        return false
    }

    /**
     * Resume a task from the retry chain of an attempt launched past an existing work dir, as written
     * before retries continued the walk: the retry of the attempt at {@code hash} was keyed
     * {@code H(hash, 2)}, then {@code H(H(hash, 2), 3)}, and so on, which the walk never reaches.
     */
    private static boolean resumeLegacyRetry(TaskRun task, HashCode hash, TaskResolver resolver) {
        for( int i = 2; i <= MAX_LOOKAHEAD + 1; i++ ) {
            hash = HashBuilder.defaultHasher().putBytes(hash.asBytes()).putInt(i).hash()
            final entry = resolver.entry(hash)
            if( !entry )
                return false
            final workDir = FileHelper.asPath(entry.trace.getWorkDir())
            if( entry.trace.isCompleted() && workDir.exists() )
                return resolver.resume(task, hash, workDir, entry)
        }
        return false
    }

}
