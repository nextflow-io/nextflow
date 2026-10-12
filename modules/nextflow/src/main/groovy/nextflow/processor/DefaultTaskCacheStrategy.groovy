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

    @Override
    boolean isEnabled(Session session) { return true }

    @Override
    void resolve(TaskRun task, HashCode hash, boolean shouldTryCache, TaskResolver resolver) {

        // only a first attempt receives the content hash -- a retry receives the previous attempt's
        // folded hash -- and clone()/makeCopy() carry it forward across the retry chain
        if( task.contentHash == null )
            task.contentHash = hash

        if( shouldTryCache && resumeFromRealizedHash(task, resolver) )
            return

        int tries = task.failCount +1
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

            // submit task for execution
            resolver.launch( task, hash, workDir )
            break
        }

    }

    /**
     * Resume a task from the realized-hash index in a single lookup, skipping the per-attempt walk,
     * so that a task which succeeded only after failed attempts is found even when the work dirs of
     * those attempts are gone or the walk would not reach it.
     *
     * @return {@code true} if the task was resumed, {@code false} to fall back to the walk (index miss,
     *      stale pointer, or lookup error). A stale pointer is rewritten by the next successful resume
     *      or execution.
     */
    private boolean resumeFromRealizedHash(TaskRun task, TaskResolver resolver) {
        try {
            final indexed = resolver.realizedHash(task.contentHash)
            if( !indexed )
                return false
            final entry = resolver.entry(indexed)
            final resumeDir = entry ? FileHelper.asPath(entry.trace.getWorkDir()) : null
            if( !entry || !entry.trace.isCompleted() || !resumeDir?.exists() )
                return false
            // the resumed copy inherits the flag, so notifyTaskCached skips rewriting a correct pointer
            task.resumedFromIndex = true
            try {
                return resolver.resume(task, indexed, resumeDir, entry)
            }
            finally {
                task.resumedFromIndex = false
            }
        }
        catch( Throwable t ) {
            log.trace "[${task.lazyName()}] Realized-hash index lookup failed -- falling back to scan -- ${t.message}"
            return false
        }
    }

}
