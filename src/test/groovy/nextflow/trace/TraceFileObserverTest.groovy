/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.trace
import java.nio.file.Files

import nextflow.executor.NopeTaskHandler
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.util.CacheHelper
import nextflow.util.Duration
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TraceFileObserverTest extends Specification {


    def testFileTrace() {

        given:
        def testFolder = Files.createTempDirectory('trace-dir')
        def file = testFolder.resolve('trace')

        // the handler
        def task = new TaskRun(id:111, name:'simple_task', hash: CacheHelper.hasher(1).hash())
        def handler = new NopeTaskHandler(task)
        def now = System.currentTimeMillis()

        // the observer class under test
        def observer = new TraceFileObserver(file)

        when:
        observer.onFlowStart(null)
        then:
        observer.current == [:]

        when:
        handler.status = TaskHandler.Status.SUBMITTED
        observer.onProcessSubmit( handler )
        def record = observer.current.get(111)
        then:
        observer.delim == '\t'
        record.taskId == 111
        record.name == 'simple_task'
        record.submit >= now
        record.start == 0
        record.complete == 0
        observer.current.containsKey(111)

        when:
        sleep 50
        handler.status = TaskHandler.Status.RUNNING
        observer.onProcessStart( handler )
        record = observer.current.get(111)
        then:
        record.start >= record.submit
        record.complete == 0
        observer.current.containsKey(111)

        when:
        sleep 50
        handler.status = TaskHandler.Status.COMPLETED
        handler.task.exitStatus = 127
        observer.onProcessComplete(handler)
        then:
        !(observer.current.containsKey(111))

        when:
        record = handler.getTraceRecord()
        observer.onFlowComplete()
        def head = file.text.readLines()[0].trim()
        def parts = file.text.readLines()[1].split('\t') as List
        then:
        head == observer.HEAD.join('\t')
        parts[0] == '111'                           // task-id
        parts[1] == 'fe/ca28af'                     // hash
        parts[2] == '-'                             // native-id
        parts[3] == 'simple_task'                   // process name
        parts[4] == TaskHandler.Status.COMPLETED.toString()
        parts[5] == '127'                           // exist-status
        TraceFileObserver.FMT.parse(parts[6]).time == record.submit         // submit time
        TraceFileObserver.FMT.parse(parts[7]).time == record.start          // start time
        TraceFileObserver.FMT.parse(parts[8]).time == record.complete       // complete time
        new Duration(parts[9]).toMillis() == record.complete -record.submit                 // wall-time
        new Duration(parts[10]).toMillis() == record.complete -record.start                 // run-time


        cleanup:
        testFolder.deleteDir()

    }

    def testBytes() {
        expect:
        TraceFileObserver.bytes(null) == '-'
        TraceFileObserver.bytes('0') == '0'
        TraceFileObserver.bytes('100') == '100 B'
        TraceFileObserver.bytes('1024') == '1 KB'
        TraceFileObserver.bytes(' 2048 ') == '2 KB'
        TraceFileObserver.bytes('abc') == 'abc'
    }

    def testTime() {
        expect:
        TraceFileObserver.time(0) == '0ms'
        TraceFileObserver.time(100) == '100ms'
        TraceFileObserver.time(2000) == '2s'
        TraceFileObserver.time(3600 * 1000 * 3) == '3h'
        TraceFileObserver.time(3600 * 1000 * 3 + 5000) == '3h 5s'
    }


    def testRender() {

        given:
        def record = new TraceRecord(
                taskId: 30,
                hash: '43d7ef',
                nativeId: '2000',
                name: 'hello',
                status: TaskHandler.Status.SUBMITTED,
                exit: 99,
                start: 1408714875000,
                submit: 1408714874000,
                complete: 1408714912000,
                cpu: '17.5033',
                mem: '3720851456.0000')
        when:
        TraceFileObserver.FMT.setTimeZone(TimeZone.getTimeZone('UTC')) // note: set the timezone to be sure the time string does not change on CI test servers
        def trace = [:] as TraceFileObserver
        def result = trace.render(record).split('\t')
        then:
        result[0] == '30'                       // task id
        result[1] == '43d7ef'                   // hash
        result[2] == '2000'                     // native id
        result[3] == 'hello'                    // name
        result[4] == 'SUBMITTED'                // status
        result[5] == '99'                       // exit status
        result[6] == '2014-08-22 13:41:14.000'  // submit
        result[7] == '2014-08-22 13:41:15.000'  // start
        result[8] == '2014-08-22 13:41:52.000'  // completed
        result[9] == '38s'                    // wall-time
        result[10] == '37s'                   // run-time
        result[11] == '17.5033'                 // cpu
        result[12] == '3.5 GB'                  // mem

    }

}
