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

package nextflow.trace
import java.nio.file.Files

import nextflow.Session
import nextflow.executor.Executor
import nextflow.executor.NopeTaskHandler
import nextflow.processor.TaskConfig
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.util.CacheHelper
import nextflow.util.Duration
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TraceFileObserverTest extends Specification {

    final static long MB = 1024 * 1024

    def setupSpec() {
        TraceRecord.TIMEZONE = TimeZone.getTimeZone('UTC') // note: set the timezone to be sure the time string does not change on CI test servers
    }

    def 'test set fields'() {

        when:
        def trace = [:] as TraceFileObserver
        trace.fields = ['task_id','name','status']
        then:
        trace.fields ==  ['task_id','name','status']

        when:
        trace = [:] as TraceFileObserver
        trace.fields = ['task_id','name','status','xxx']
        then:
        thrown(IllegalArgumentException)

    }


    def 'test set formats'() {

        when:
        def trace = [:] as TraceFileObserver
        trace.formats = ['str','num','date']
        then:
        trace.formats == ['str','num','date']

    }

    def 'test set fields and formats'() {

        def trace

        when:
        trace = [:] as TraceFileObserver
        trace.setFieldsAndFormats(['task_id:str','name:str','status:num', 'start', 'duration'])
        then:
        trace.fields ==  ['task_id','name','status', 'start', 'duration']
        trace.formats == ['str','str','num', 'date', 'time']

        when:
        trace = [:] as TraceFileObserver
        trace.useRawNumbers(true)
        trace.setFieldsAndFormats(['task_id:str','name:str','status:num', 'start:date', 'duration:time', 'realtime'])
        then:
        trace.fields ==  ['task_id','name','status', 'start', 'duration', 'realtime']
        trace.formats == ['str','str','num', 'num', 'num', 'num']

    }

    def 'test file trace'() {

        given:
        def testFolder = Files.createTempDirectory('trace-dir')
        def file = testFolder.resolve('trace')

        // the handler
        def task = new TaskRun(id:TaskId.of(111), name:'simple_task', hash: CacheHelper.hasher(1).hash(), config: new TaskConfig())
        task.processor = Mock(TaskProcessor)
        task.processor.getSession() >> new Session()
        task.processor.getName() >> 'x'
        task.processor.getExecutor() >> Mock(Executor)
        task.processor.getProcessEnvironment() >> [:]

        def handler = new NopeTaskHandler(task)
        def now = System.currentTimeMillis()

        // the observer class under test
        def observer = new TraceFileObserver(file)

        when:
        observer.onFlowCreate(null)
        then:
        observer.current.isEmpty()

        when:
        handler.status = TaskStatus.SUBMITTED
        observer.onProcessSubmit( handler, handler.getTraceRecord() )
        def record = observer.current.get(TaskId.of(111))
        then:
        observer.separator == '\t'
        record.taskId == 111L
        record.name == 'simple_task'
        record.submit >= now
        record.start == 0
        observer.current.containsKey(TaskId.of(111))

        when:
        sleep 50
        handler.status = TaskStatus.RUNNING
        observer.onProcessStart( handler, handler.getTraceRecord() )
        record = observer.current.get(TaskId.of(111))
        then:
        record.start >= record.submit
        observer.current.containsKey(TaskId.of(111))

        when:
        sleep 50
        handler.status = TaskStatus.COMPLETED
        handler.task.exitStatus = 127
        observer.onProcessComplete(handler, handler.getTraceRecord())
        then:
        !(observer.current.containsKey(TaskId.of(111)))

        when:
        record = handler.getTraceRecord()
        observer.onFlowComplete()
        def head = file.text.readLines()[0].trim()
        def parts = file.text.readLines()[1].split('\t') as List
        then:
        head == observer.fields.join('\t')
        parts[0] == '111'                           // task-id
        parts[1] == 'fe/ca28af'                     // hash
        parts[2] == '-'                             // native-id
        parts[3] == 'simple_task'                   // process name
        parts[4] == TaskStatus.COMPLETED.toString()
        parts[5] == '127'                           // exist-status
        parts[6] == TraceRecord.fmtDate(record.submit,null) // submit time
        parts[7] == new Duration(record.complete -record.submit).toString()         // wall-time
        parts[8] == new Duration(record.complete -record.start).toString()         // run-time

        cleanup:
        testFolder.deleteDir()

    }


    def 'test render'() {

        given:
        def record = new TraceRecord()
        record.task_id = 30
        record.hash = '43d7ef'
        record.native_id = '2000'
        record.name = 'hello'
        record.status = TaskStatus.COMPLETED
        record.exit = 99
        record.start = 1408714875000
        record.submit = 1408714874000
        record.complete = 1408714912000
        record.duration = 1408714912000 - 1408714874000
        record.realtime = 1408714912000 - 1408714875000
        record.'%cpu' = 17.50f
        record.peak_rss = 10_000 * 1024
        record.peak_vmem = 30_000 * 1024
        record.rchar = 30_000 * 1024
        record.wchar = 10_000 * 1024

        when:
        def trace = [:] as TraceFileObserver
        def result = trace.render(record).split('\t')
        then:
        result[0] == '30'                       // task id
        result[1] == '43d7ef'                   // hash
        result[2] == '2000'                     // native id
        result[3] == 'hello'                    // name
        result[4] == 'COMPLETED'                // status
        result[5] == '99'                       // exit status
        result[6] == '2014-08-22 13:41:14.000'  // submit
        result[7] == '38s'                      // wall-time
        result[8] == '37s'                      // run-time
        result[9] == '17.5%'                    // cpu
        result[10] == '9.8 MB'                  // peak_rss
        result[11] == '29.3 MB'                 // peak_vmem
        result[12] == '29.3 MB'                 // rchar
        result[13] == '9.8 MB'                  // wchar

    }

    def 'test custom render' () {


        given:
        def record = new TraceRecord()
        record.task_id = '5'
        record.syscr = 10
        record.syscw = 20
        record.rss = 10 * MB

        when:
        def trace = [:] as TraceFileObserver
        trace.setFieldsAndFormats( 'task_id,syscr,syscw,rss,rss:num' )
        def result = trace.render(record).split('\t')
        then:
        result[0] == '5'
        result[1] == '10'
        result[2] == '20'
        result[3] == '10 MB'
        result[4] == '10485760'

    }

    def 'should create a record in the trace file'() {

        given:
        final KB = 1024L
        final file = TestHelper.createInMemTempFile('trace')
        file.text =  '''
                pid state %cpu %mem vmem rss peak_vmem peak_rss rchar wchar syscr syscw read_bytes write_bytes
                18 0 7999 46 7868980 6694900 7876620 6702144 2147483647 2147483647 44001533 148401890 2147483647 2147483647
                9005022
                '''
                .leftTrim()

        when:
        def handler = [:] as TraceRecord
        def record = handler.parseTraceFile(file)
        record.task_id = 3
        record.hash = '9a/a894b2'
        record.native_id = '2000'
        record.name = 'mapping'
        record.status = TaskStatus.COMPLETED
        record.exit = 0
        record.start = 1429031453000
        record.submit =1429031454000
        record.complete = 1429031455000
        record.duration = 9087193
        record.realtime = 9005022
        record.queue = 'bigjobs'

        def trace = [:] as TraceFileObserver
        trace.setFieldsAndFormats('task_id,hash,native_id,name,status,exit,submit,duration,realtime,%cpu,rss,vmem,peak_rss,peak_vmem,rchar,wchar,syscr,syscw,duration:num,realtime:num,rss:num,vmem:num,peak_rss:num,peak_vmem:num,rchar:num,wchar:num,queue')
        def result = trace.render(record).split('\t')

        then:
        result[0] == '3'
        result[1] == '9a/a894b2'
        result[2] == '2000'
        result[3] == 'mapping'
        result[4] == 'COMPLETED'        //status
        result[5] == '0'                // exit
        result[6] == '2015-04-14 17:10:54.000'    // submit
        result[7] == '2h 31m 27s'       // duration
        result[8] == '2h 30m 5s'        // realtime
        result[9] == '799.9%'           // %cpu
        result[10] == '6.4 GB'          // rss
        result[11] == '7.5 GB'          // vmem
        result[12] == '6.4 GB'          // peak_rss
        result[13] == '7.5 GB'          // peak_vmem
        result[14] == '2 GB'            // rchar
        result[15] == '2 GB'            // wchar
        result[16] == '44001533'        // syscr
        result[17] == '148401890'       // syscw
        result[18] == '9087193'         // duration:num
        result[19] == '9005022'         // realtime:num
        result[20] == 6694900*KB as String      // rss:num
        result[21] == 7868980*KB as String      // vmem:num
        result[22] == 6702144*KB as String      // peak_rss:num
        result[23] == 7876620*KB as String      // peak_vmem:num
        result[24] == '2147483647'    // rchar:num
        result[25] == '2147483647'    // wchar:num
        result[26] == 'bigjobs'       // queue
    }

}
