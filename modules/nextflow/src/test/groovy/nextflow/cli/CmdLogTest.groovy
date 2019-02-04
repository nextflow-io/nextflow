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

package nextflow.cli
import java.nio.file.Files

import nextflow.CacheDB
import nextflow.executor.CachedTaskHandler
import nextflow.processor.ProcessConfig
import nextflow.processor.TaskContext
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.TaskBody
import nextflow.trace.TraceRecord
import nextflow.util.CacheHelper
import nextflow.util.HistoryFile
import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdLogTest extends Specification {

    /*
     * Read more http://mrhaki.blogspot.com.es/2015/02/spocklight-capture-and-assert-system.html
     */
    @Rule
    OutputCapture capture = new OutputCapture()

    def 'should print the folders for executed tasks' () {

        setup:
        final folder = Files.createTempDirectory('test')
        folder.resolve('.nextflow').mkdir()
        final uuid = UUID.randomUUID()
        final runName = 'test_1'

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
        def task1 = Mock(TaskRun)
        task1.getProcessor() >> proc
        task1.getHash() >> { CacheHelper.hasher('x1').hash() }

        def task2 = Mock(TaskRun)
        task2.getProcessor() >> proc
        task2.getHash() >> { CacheHelper.hasher('x2').hash() }

        def task3 = Mock(TaskRun)
        task3.getProcessor() >> proc
        task3.getHash() >> { CacheHelper.hasher('x3').hash() }

        cache.open()
        def h1 = new CachedTaskHandler(task1, new TraceRecord([task_id: 1, process: 'foo', exit: 0, workdir:"$folder/aaa"]))
        cache.writeTaskIndex0(h1)
        cache.writeTaskEntry0(h1, h1.getTraceRecord())

        def h2 = new CachedTaskHandler(task2, new TraceRecord([task_id: 2, process: 'foo', exit: 1, workdir:"$folder/bbb"]))
        cache.writeTaskIndex0(h2)
        cache.writeTaskEntry0(h2, h2.getTraceRecord())

        def h3 = new CachedTaskHandler(task3, new TraceRecord([task_id: 3, process: 'bar', exit: 0, workdir:"$folder/ccc"]))
        cache.writeTaskIndex0(h3)
        cache.writeTaskEntry0(h3, h3.getTraceRecord())
        cache.close()

        def history = new HistoryFile(folder.resolve(HistoryFile.FILE_NAME))
        history.write(runName,uuid,'b3d3aca8eb','run')

        when:
        def log = new CmdLog(basePath: folder, args: [runName])
        log.run()
        def stdout = capture
                .toString()
                .readLines()
                // remove the log part
                .findResults { line -> !line.contains('DEBUG') ? line : null }
                .join('\n')
        then:
        stdout.readLines().size() == 3
        stdout.readLines().contains( "$folder/aaa" .toString())
        stdout.readLines().contains( "$folder/bbb" .toString())
        stdout.readLines().contains( "$folder/ccc" .toString())

    }

    def 'should filter printed tasks' () {

        setup:
        final folder = Files.createTempDirectory('test')
        folder.resolve('.nextflow').mkdir()
        final uuid = UUID.randomUUID()
        final runName = 'test_1'

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
        def task1 = Mock(TaskRun)
        task1.getProcessor() >> proc
        task1.getHash() >> { CacheHelper.hasher('x1').hash() }

        def task2 = Mock(TaskRun)
        task2.getProcessor() >> proc
        task2.getHash() >> { CacheHelper.hasher('x2').hash() }

        def task3 = Mock(TaskRun)
        task3.getProcessor() >> proc
        task3.getHash() >> { CacheHelper.hasher('x3').hash() }

        cache.open()
        def h1 = new CachedTaskHandler(task1, new TraceRecord([task_id: 1, process: 'foo', exit: 0, workdir:"$folder/aaa"]))
        cache.writeTaskIndex0(h1)
        cache.writeTaskEntry0(h1, h1.getTraceRecord())

        def h2 = new CachedTaskHandler(task2, new TraceRecord([task_id: 2, process: 'foo', exit: 1, workdir:"$folder/bbb"]))
        cache.writeTaskIndex0(h2)
        cache.writeTaskEntry0(h2, h2.getTraceRecord())

        def h3 = new CachedTaskHandler(task3, new TraceRecord([task_id: 3, process: 'bar', exit: 0, workdir:"$folder/ccc"]))
        cache.writeTaskIndex0(h3)
        cache.writeTaskEntry0(h3, h3.getTraceRecord())
        cache.close()

        def history = new HistoryFile(folder.resolve(HistoryFile.FILE_NAME))
        history.write(runName,uuid,'b3d3aca8eb','run')


        when:
        def log = new CmdLog(basePath: folder, filterStr: 'exit == 0', args: ['test_1'])
        log.run()

        def stdout = capture
                .toString()
                .readLines()
                // remove the log part
                .findResults { line -> !line.contains('DEBUG') ? line : null }
                .join('\n')
        then:
        stdout.readLines().size() == 2
        stdout.readLines().contains( "$folder/aaa" .toString())
        stdout.readLines().contains( "$folder/ccc" .toString())

    }

}
