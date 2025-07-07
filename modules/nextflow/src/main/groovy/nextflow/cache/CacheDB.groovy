/*
 * Copyright 2013-2024, Seqera Labs
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
 *
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

    void incTaskEntry( HashCode hash ) {
        final payload = store.getEntry(hash)
        if( !payload ) {
            log.debug "Can't increment reference for cached task with key: $hash"
            return
        }

        final record = (List)KryoHelper.deserialize(payload)
        // third record contains the reference count for this record
        record[2] = ((Integer)record[2]) +1
        // save it again
        store.putEntry(hash, KryoHelper.serialize(record))

    }

    boolean removeTaskEntry( HashCode hash ) {
        final payload = store.getEntry(hash)
        if( !payload ) {
            log.debug "Can't increment reference for cached task with key: $hash"
            return false
        }

        final record = (List)KryoHelper.deserialize(payload)
        // third record contains the reference count for this record
        def count = record[2] = ((Integer)record[2]) -1
        // save or delete
        if( count > 0 ) {
            store.putEntry(hash, KryoHelper.serialize(record))
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
        writer.send { writeTaskEntry0(handler, trace) }
    }

    void cacheTaskAsync( TaskHandler handler ) {
        writer.send {
            writeTaskIndex0(handler,true)
            incTaskEntry(handler.task.hash)
        }
    }

    void putIndexAsync(TaskHandler handler ) {
        writer.send { writeTaskIndex0(handler) }
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
