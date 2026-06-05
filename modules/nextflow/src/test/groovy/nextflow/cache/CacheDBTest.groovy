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

import com.google.common.hash.HashCode
import nextflow.executor.CachedTaskHandler
import nextflow.processor.TaskContext
import nextflow.processor.TaskEntry
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import nextflow.trace.TraceRecord
import nextflow.util.CacheHelper
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

    def 'should delegate hash index read and write' () {
        given:
        def store = Mock(CacheStore)
        def cache = new CacheDB(store)
        and:
        def content = HashCode.fromInt(1)
        def finalHash = HashCode.fromInt(2)

        when:
        def result = cache.getSuccessfulHash(content)
        then:
        1 * store.getSuccessfulHash(content) >> finalHash
        result == finalHash

        when:
        cache.putSuccessfulHashAsync(content, finalHash)
        cache.close()
        then:
        1 * store.putSuccessfulHash(content, finalHash)
    }

    private makeHandler(HashCode hash, HashCode contentHash) {
        def proc = Mock(TaskProcessor)
        proc.getTaskBody() >> new BodyDef(null,'source')
        proc.getConfig() >> new ProcessConfig([:])
        def task = Mock(TaskRun)
        task.getProcessor() >> proc
        task.getHash() >> hash
        task.getContentHash() >> contentHash
        return new CachedTaskHandler(task, new TraceRecord())
    }

    def 'should remove the hash index pointer when its target entry is deleted' () {
        given:
        def folder = Files.createTempDirectory('test')
        def store = new DefaultCacheStore(UUID.randomUUID(), 'test_1', folder)
        def cache = new CacheDB(store).open()
        and:
        def content = CacheHelper.hasher('CONTENT').hash()
        def successHash = CacheHelper.hasher('SUCCESS').hash()
        def failHash = CacheHelper.hasher('FAIL').hash()
        and: 'a successful entry carrying its content hash, indexed'
        cache.writeTaskEntry0(makeHandler(successHash, content), new TraceRecord([task_id:1, exit:0]))
        store.putSuccessfulHash(content, successHash)
        and: 'a failed-attempt entry for the SAME content hash (not the index target)'
        cache.writeTaskEntry0(makeHandler(failHash, content), new TraceRecord([task_id:1, exit:1]))

        when: 'the failed-attempt entry is removed'
        def removedFail = cache.removeTaskEntry(failHash)
        then: 'the entry is gone but the pointer to the successful entry survives'
        removedFail
        store.getSuccessfulHash(content) == successHash

        when: 'the successful entry (the index target) is removed'
        def removedSuccess = cache.removeTaskEntry(successHash)
        then: 'its index pointer is removed too'
        removedSuccess
        store.getSuccessfulHash(content) == null

        cleanup:
        cache?.close()
        folder?.deleteDir()
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

}
