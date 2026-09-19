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

import java.nio.file.Files

import ch.qos.logback.classic.Level
import ch.qos.logback.classic.Logger
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.core.read.ListAppender
import com.google.common.hash.HashCode
import nextflow.cache.CacheDB
import nextflow.cache.DefaultCacheStore
import nextflow.executor.CachedTaskHandler
import nextflow.processor.TaskContext
import nextflow.processor.TaskEntry
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import nextflow.trace.TraceRecord
import nextflow.util.CacheHelper
import nextflow.util.KryoHelper
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CacheDBTest extends Specification {


    def 'should save and read a task entry in the cache db' () {

        setup:
        def folder = Files.createTempDirectory('test')
        def uuid = UUID.randomUUID()
        def hash = CacheHelper.hasher('x').hash()
        def runName = 'test_1'

        // -- the session object
        def store = new DefaultCacheStore(uuid, runName, folder)
        def cache = new CacheDB(store)

        // -- the processor mock
        def proc = Mock(TaskProcessor)
        proc.getTaskBody() >> new BodyDef(null,'source')
        proc.getConfig() >> new ProcessConfig([:])

        // -- the task context
        def ctx = new TaskContext()
        ctx.setHolder( [X: 10, Y: 'Hello'] )

        // -- the task mock
        def task = Mock(TaskRun)
        task.getProcessor() >> proc
        task.getHash() >> hash
        task.getId() >> TaskId.of(2)

        when:
        cache.open()
        then:
        folder.resolve("cache/$uuid/db").exists()
        folder.resolve("cache/$uuid/index.$runName").exists()

        when:
        def trace = new TraceRecord([task_id: 1, process: 'foo', exit: 0])
        def handler = new CachedTaskHandler(task, trace)
        cache.writeTaskEntry0( handler, trace )
        then:
        1 * proc.isCacheable() >> true
        1 * task.hasCacheableValues() >> true
        1 * task.getContext() >> ctx

        when:
        def entry = cache.getTaskEntry(hash, proc)
        then:
        entry instanceof TaskEntry
        entry.trace instanceof TraceRecord
        entry.trace.get('task_id') == 2   // task_id is taken from the current TaskRun
        entry.trace.get('process') == 'foo'
        entry.trace.get('exit') == 0
        entry.context instanceof TaskContext
        entry.context.X == 10
        entry.context.Y == 'Hello'

        cleanup:
        cache?.close()
        folder?.deleteDir()
    }


    private makeTaskHandler(HashCode hash, Map record, Map context=null) {

        // -- the processor mock
        def proc = Mock(TaskProcessor)
        proc.getTaskBody() >> new BodyDef(null,'source')
        proc.getConfig() >> new ProcessConfig([:])

        // -- the task context
        def ctx = new TaskContext()
        if( context )
            ctx.setHolder(context)

        // -- the task mock
        def task = Mock(TaskRun)
        task.getProcessor() >> proc
        task.getHash() >> hash

        def trace = new TraceRecord()
        return new CachedTaskHandler(task, trace)

    }

    def 'should write some tasks and iterate over them' () {

        setup:
        def folder = Files.createTempDirectory('test')
        def uuid = UUID.randomUUID()
        def hash1 = CacheHelper.hasher('x').hash()
        def hash2 = CacheHelper.hasher('x').hash()
        def hash3 = CacheHelper.hasher('x').hash()
        def runName = 'test_1'
        def trace = Mock(TraceRecord)
        and:
        def store = new DefaultCacheStore(uuid, runName, folder)

        when:
        def cache = new CacheDB(store).open()

        def h1 = makeTaskHandler(hash1, [task_id: 1, process: 'foo', exit: 0])
        cache.writeTaskEntry0(h1, h1.traceRecord)
        cache.writeTaskIndex0(h1)

        def h2 = makeTaskHandler(hash2, [task_id: 2, process: 'bar', exit: 0])
        cache.writeTaskEntry0(h2, h1.traceRecord)
        cache.writeTaskIndex0(h2)

        def h3 = makeTaskHandler(hash3, [task_id: 3, process: 'baz', exit: 1])
        cache.writeTaskEntry0(h3, h1.traceRecord)
        cache.writeTaskIndex0(h3)

        // done
        cache.close()

        then:
        noExceptionThrown()


        when:
        cache.openForRead()
        def items = []
        cache.eachRecord { k, v -> items << ( [hash: k, record: v] ) }
        cache.close()
        then:
        items.size() == 3
        items[0].hash == hash1
        items[1].hash == hash2
        items[2].hash == hash3


        cleanup:
        folder?.deleteDir()

    }

    def 'the refcount read-modify-write goes through updateEntry, a new entry through putEntry' () {
        given:
        def hash = CacheHelper.hasher('x').hash()
        def record = [new TraceRecord([task_id:1]).serialize(), null, 1]
        def store = Mock(CacheStore)
        def cache = new CacheDB(store)

        when: 'the reference count of an existing entry is bumped'
        cache.incTaskEntry(hash)
        then: 'it is an update of the record just read -- a composite store must not treat it as new'
        1 * store.getEntry(hash) >> KryoHelper.serialize(record)
        1 * store.updateEntry(hash, _)
        0 * store.putEntry(_, _)

        when: 'the reference count is decremented but the entry survives'
        cache.removeTaskEntry(hash)
        then:
        1 * store.getEntry(hash) >> KryoHelper.serialize([record[0], null, 2])
        1 * store.updateEntry(hash, _)
        0 * store.putEntry(_, _)

        when: 'a brand new entry is recorded'
        def proc = Mock(TaskProcessor) { getConfig() >> new ProcessConfig([:]) }
        def task = Mock(TaskRun) { getProcessor() >> proc; getHash() >> hash }
        cache.writeTaskEntry0(new CachedTaskHandler(task, new TraceRecord()), new TraceRecord([task_id:1]))
        then: 'it must go to the writable store, so it is a plain putEntry'
        1 * store.putEntry(hash, _)
        0 * store.updateEntry(_, _)
    }

    def 'updateEntry defaults to putEntry for a store that does not override it' () {
        given:
        def store = Spy(DefaultCacheStore, constructorArgs: [UUID.randomUUID(), 'r', Files.createTempDirectory('test')])
        def hash = CacheHelper.hasher('x').hash()

        when:
        store.open()
        store.updateEntry(hash, 'hello'.bytes)
        then:
        1 * store.putEntry(hash, _)
        and:
        new String(store.getEntry(hash)) == 'hello'

        cleanup:
        store?.close()
    }


    def 'a failed async cache write is logged, not silently swallowed' () {
        given: 'a store whose write throws (e.g. a transient cloud error)'
        def store = Stub(CacheStore) {
            writeIndex(_, _) >> { throw new RuntimeException('s3 boom') }
        }
        def cache = new CacheDB(store)
        and: 'capture CacheDB logs'
        def logger = (Logger) org.slf4j.LoggerFactory.getLogger(CacheDB)
        def appender = new ListAppender<ILoggingEvent>()
        appender.start()
        logger.addAppender(appender)

        when:
        cache.putIndexAsync(Mock(TaskHandler) { getTask() >> new TaskRun(hash: HashCode.fromInt(1)) })
        cache.close()   // close() awaits the writer agent, so the failing write has already run

        then: 'the failure surfaced as a WARN instead of being swallowed by the agent'
        appender.list.any { it.level == Level.WARN && it.formattedMessage.contains('Unable to persist cache record') }

        cleanup:
        logger.detachAppender(appender)
    }

}
