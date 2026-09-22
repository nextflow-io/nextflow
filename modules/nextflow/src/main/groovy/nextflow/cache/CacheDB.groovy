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

package nextflow.cache


import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.SimpleType
import groovy.util.logging.Slf4j
import groovyx.gpars.agent.Agent
import nextflow.processor.TaskContext
import nextflow.processor.TaskEntry
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
import nextflow.trace.TraceRecord
import nextflow.util.KryoHelper
/**
 * Manages nextflow cache DB
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CacheDB implements Closeable {

    /** An agent used to apply asynchronously DB write operations */
    private Agent writer

    private CacheStore store

    CacheDB(CacheStore store) {
        this.store = store
        this.writer = new Agent()
    }

    /**
     * Initialise the database structure on the underlying file system
     *
     * @return The {@link CacheDB} instance itself
     */
    CacheDB open() {
        store.open()
        return this
    }

    /**
     * Open the database in read mode
     *
     * @return The {@link CacheDB} instance itself
     */
    CacheDB openForRead() {
        store.openForRead()
        return this
    }

    /**
     * Retrieve a task runtime information from the cache DB
     *
     * @param taskHash The {@link HashCode} of the task to retrieve
     * @param processor The {@link TaskProcessor} instance to be assigned to the retrieved task
     * @return A {link TaskEntry} instance or {@code null} if a task for the given hash does not exist
     */
    TaskEntry getTaskEntry(HashCode taskHash, TaskProcessor processor) {

        final payload = store.getEntry(taskHash)
        if( !payload )
            return null

        final record = (List)KryoHelper.deserialize(payload)
        TraceRecord trace = TraceRecord.deserialize( (byte[])record[0] )
        TaskContext ctx = record[1]!=null && processor!=null ? TaskContext.deserialize(processor, (byte[])record[1]) : null

        return new TaskEntry(trace,ctx)
    }

    /**
     * Dispatch an asynchronous cache write on the {@link #writer} agent, logging any failure instead
     * of letting the agent swallow it. Without this the run still succeeds while the cache is
     * silently incomplete, with no trace of why.
     *
     * <p>What a dropped write costs depends on which one it was, so the message stays neutral:
     * a failed {@code putTaskAsync} leaves the task non-resumable, so it re-executes on the next
     * resume; a failed {@code putIndexAsync} only affects {@code nextflow log} and {@code clean},
     * since resume reads the entry through {@link #getTaskEntry}; a failed {@code cacheTaskAsync}
     * leaves the reference count stale. Surfacing the error makes all three visible (and is the
     * natural hook for a future retry at the store level).
     *
     * @param what   short description of the write for the log message (typically the task hash)
     * @param action the store mutation to run on the writer thread
     */
    protected void dispatchWrite(String what, Closure action) {
        writer.send {
            try {
                action.call()
            }
            catch( Throwable e ) {
                log.warn("Unable to persist cache record for ${what}", e)
            }
        }
    }

    /**
     * Bump the reference count of an entry that already exists, which also refreshes the stored
     * object's last-modified stamp.
     *
     * <b>The count is advisory for a shared store.</b> This is a read-modify-write with no
     * compare-and-set behind {@link CacheStore}, so two runs resuming the same task can read the same
     * value and write the same increment, losing one. That is harmless today only because neither
     * consumer of the count is reachable for a shared cache: {@link #removeTaskEntry} — the only
     * decrement, and the only path that can delete on reaching zero — is called from
     * {@code Session.cleanup}, which returns before opening the cache for a non-{@code file:} work
     * dir, and from {@code CmdClean}, which a shared cache's {@code CacheDB} may refuse by
     * overriding it; and a cross-run cache ages entries by the object's last-modified time, not by
     * the count. Eviction tooling that reads the count instead would need a genuinely atomic update
     * here.
     */
    void incTaskEntry( HashCode hash ) {
        final payload = store.getEntry(hash)
        if( !payload ) {
            log.debug "Can't increment reference for cached task with key: $hash"
            return
        }

        final record = (List)KryoHelper.deserialize(payload)
        // third record contains the reference count for this record
        record[2] = ((Integer)record[2]) +1
        // save it again -- an update of the record just read, not a new entry (see updateEntry)
        store.updateEntry(hash, KryoHelper.serialize(record))

    }

    /**
     * Decrement the reference count, deleting the entry when it reaches zero. Callers must not invoke
     * this on a shared cache — see the advisory-count note on {@link #incTaskEntry}.
     */
    boolean removeTaskEntry( HashCode hash ) {
        final payload = store.getEntry(hash)
        if( !payload ) {
            log.debug "Can't increment reference for cached task with key: $hash"
            return false
        }

        final record = (List)KryoHelper.deserialize(payload)
        // third record contains the reference count for this record
        def count = record[2] = ((Integer)record[2]) -1
        // save or delete -- as in incTaskEntry, saving is an update of the record just read
        if( count > 0 ) {
            store.updateEntry(hash, KryoHelper.serialize(record))
            return false
        }
        else {
            store.deleteEntry(hash)
            return true
        }
    }


    /**
     * Save task runtime information to th cache DB
     *
     * @param handler A {@link TaskHandler} instance
     */
    @PackageScope
    void writeTaskEntry0( TaskHandler handler, TraceRecord trace ) {

        final task = handler.task
        final proc = task.processor
        final key = task.hash

        // save the context map for caching purpose
        // only the 'cache' is active and
        TaskContext ctx = proc.isCacheable() && task.hasCacheableValues() ? task.context : null

        def record = new ArrayList(3)
        record[0] = trace.serialize()
        record[1] = ctx != null ? ctx.serialize() : null
        record[2] = 1

        // -- save in the db
        store.putEntry( key, KryoHelper.serialize(record) )

    }

    void putTaskAsync( TaskHandler handler, TraceRecord trace ) {
        dispatchWrite("task entry ${handler.task.hash}") { writeTaskEntry0(handler, trace) }
    }

    void cacheTaskAsync( TaskHandler handler ) {
        dispatchWrite("cached task ${handler.task.hash}") {
            writeTaskIndex0(handler,true)
            incTaskEntry(handler.task.hash)
        }
    }

    void putIndexAsync(TaskHandler handler ) {
        dispatchWrite("task index ${handler.task.hash}") { writeTaskIndex0(handler) }
    }

    @PackageScope
    void writeTaskIndex0( TaskHandler handler, boolean cached = false ) {
        store.writeIndex(handler.task.hash, cached)
    }

    void deleteIndex() {
        store.deleteIndex()
    }

    void drop() {
        store.drop()
    }

    /**
     * Iterate the tasks cache using the index file
     * @param closure The operation to applied
     * @return The {@link CacheDB} instance itself
     */
    CacheDB eachRecord( Closure closure ) {
        assert closure

        final itr = store.iterateIndex()
        while( itr.hasNext() ) {
            final index = itr.next()

            final payload = store.getEntry(index.key)
            if( !payload ) {
                log.trace "Unable to retrieve cache record for key: ${-> index.key}"
                continue
            }

            final record = (List<byte[]>)KryoHelper.deserialize(payload)
            TraceRecord trace = TraceRecord.deserialize(record[0])
            trace.setCached(index.cached)

            final refCount = record[2] as Integer

            final len=closure.maximumNumberOfParameters
            if( len==1 )
                closure.call(trace)

            else if( len==2 )
                closure.call(index.key, trace)

            else if( len==3 )
                closure.call(index.key, trace, refCount)

            else
                throw new IllegalArgumentException("Invalid closure signature -- Too many parameters")

        }

        return this
    }

    TraceRecord getTraceRecord( HashCode hashCode ) {
        final result = getTaskEntry(hashCode, null)
        return result ? result.trace : null
    }

    TraceRecord findTraceRecord( @ClosureParams(value = SimpleType.class, options = "nextflow.trace.TraceRecord") Closure<Boolean> criteria ) {

        final itr = store.iterateIndex()
        while( itr.hasNext() ) {
            final index = itr.next()

            final payload = store.getEntry(index.key)
            if( !payload ) {
                log.trace "Unable to retrieve cache record for key: ${-> index.key}"
                continue
            }

            final record = (List<byte[]>)KryoHelper.deserialize(payload)
            TraceRecord trace = TraceRecord.deserialize(record[0])
            trace.setCached(index.cached)

            final len=criteria.maximumNumberOfParameters
            if( len!=1 ) {
                throw new IllegalArgumentException("Invalid criteria signature -- Too many parameters")
            }
            if( criteria.call(trace) )
                return trace
        }
        // no matches
        return null
    }

    /**
     * Close the underlying database and index file
     */
    @Override
    void close() {
        log.trace "Closing CacheDB.."
        writer.await()
        log.trace "Closing CacheDB index"
        store.close()
        log.debug "Closing CacheDB done"
    }
}
