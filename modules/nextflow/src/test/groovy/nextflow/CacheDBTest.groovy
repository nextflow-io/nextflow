/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow
import java.nio.file.Files

import com.google.common.hash.HashCode
import nextflow.executor.CachedTaskHandler
import nextflow.processor.TaskId
import nextflow.script.ProcessConfig
import nextflow.processor.TaskContext
import nextflow.processor.TaskEntry
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.BodyDef
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
        def cache = new CacheDB(uuid, runName, folder)

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

        when:
        def cache = new CacheDB(uuid, runName, folder).open()

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
