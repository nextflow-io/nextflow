/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow
import java.nio.file.Files

import com.google.common.hash.HashCode
import nextflow.executor.CachedTaskHandler
import nextflow.processor.ProcessConfig
import nextflow.processor.TaskContext
import nextflow.processor.TaskEntry
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.TaskBody
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
        proc.getTaskBody() >> new TaskBody(null,'source')
        proc.getConfig() >> new ProcessConfig([:])

        // -- the task context
        def ctx = new TaskContext()
        ctx.setHolder( [X: 10, Y: 'Hello'] )

        // -- the task mock
        def task = Mock(TaskRun)
        task.getProcessor() >> proc
        task.getHash() >> hash

        when:
        cache.open()
        then:
        folder.resolve("cache/$uuid/db").exists()
        folder.resolve("cache/$uuid/index.$runName").exists()

        when:
        def trace = new TraceRecord([task_id: 1, process: 'foo', exit: 0])
        def handler = new CachedTaskHandler(task, trace)
        cache.writeTaskEntry0( handler )
        then:
        1 * proc.isCacheable() >> true
        1 * task.hasCacheableValues() >> true
        1 * task.getContext() >> ctx

        when:
        def entry = cache.getTaskEntry(hash, proc)
        then:
        entry instanceof TaskEntry
        entry.trace instanceof TraceRecord
        entry.trace.get('task_id') == 1
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
        proc.getTaskBody() >> new TaskBody(null,'source')
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

        when:
        def cache = new CacheDB(uuid, runName, folder).open()

        def h1 = makeTaskHandler(hash1, [task_id: 1, process: 'foo', exit: 0])
        cache.writeTaskEntry0(h1)
        cache.writeTaskIndex0(h1)

        def h2 = makeTaskHandler(hash2, [task_id: 2, process: 'bar', exit: 0])
        cache.writeTaskEntry0(h2)
        cache.writeTaskIndex0(h2)

        def h3 = makeTaskHandler(hash3, [task_id: 3, process: 'baz', exit: 1])
        cache.writeTaskEntry0(h3)
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
