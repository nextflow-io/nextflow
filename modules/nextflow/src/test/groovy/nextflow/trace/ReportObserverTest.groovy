/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
import java.nio.file.Paths
import java.util.concurrent.Executors

import groovy.json.JsonSlurper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskId
import nextflow.script.WorkflowMetadata
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ReportObserverTest extends Specification {

    def setupSpec() {
        TraceRecord.TIMEZONE = TimeZone.getTimeZone('UTC') // note: set the timezone to be sure the time string does not change on CI test servers
    }

    def 'should render json data' () {

        given:
        def now = 1429821425141
        def r1 = new TraceRecord()
        r1.task_id = '1'
        r1.name = 'foo'
        r1.process = 'alpha'

        def r2 = new TraceRecord()
        r2.task_id = '2'
        r2.name = 'bar'
        r2.submit = now
        r2.start = now + 100
        r2.complete = now + 500
        r2.realtime = 400
        r2.duration = 500
        r2.process = 'alpha'

        def r3 = new TraceRecord()
        r3.task_id = '3'
        r3.name = 'baz'
        r3.submit = now
        r3.start = now + 200
        r3.complete = now + 700
        r3.realtime = 500
        r3.duration = 700
        r3.process = 'beta'

        def observer = [:] as ReportObserver
        def data = [r1, r2, r3]

        when:
        def str = observer.renderJsonData(data)
        def trace = new JsonSlurper().parseText(str) as List

        then:
        trace[0].task_id == '1'
        trace[0].name == 'foo'
        trace[0].process == 'alpha'

        trace[1].task_id == '2'
        trace[1].name == 'bar'
        trace[1].submit == '1429821425141'
        trace[1].start == '1429821425241'
        trace[1].complete == '1429821425641'
        trace[1].realtime == '400'
        trace[1].duration == '500'
        trace[1].process == 'alpha'

        trace[2].task_id == '3'
        trace[2].name == 'baz'
        trace[2].submit == '1429821425141'
        trace[2].start == '1429821425341'
        trace[2].complete == '1429821425841'
        trace[2].realtime == '500'
        trace[2].duration == '700'
        trace[2].process == 'beta'

    }

    def 'should render json payload' () {
        given:
        def report = Spy(ReportObserver)

        when:
        def result = report.renderPayloadJson()
        then:
        1 * report.renderTasksJson() >> '[1,2,3]'
        1 * report.renderSummaryJson() >> '{foo: 1}'

        result == '{ "trace":[1,2,3], "summary":{foo: 1} }'
    }

    def 'should render html template' () {

        given:
        def workDir = TestHelper.createInMemTempDir()
        def workflow = new WorkflowMetadata(
                workDir: workDir,
                stats: new WorkflowStats(),
                nextflow: [version: '0.27.9', build: '3232', timestamp: '2017-12-12']
        )

        def file = TestHelper.createInMemTempFile('report.html')
        def observer = Spy(ReportObserver, constructorArgs: [file])

        when:
        observer.renderHtml()
        then:
        1 * observer.getWorkflowMetadata() >> workflow
        1 * observer.renderSummaryJson() >> '{ }'
        1 * observer.renderTasksJson() >> '{ }'
        Files.exists(file)
        file.text.contains('Nextflow Report')
        file.text.contains( workDir.toUriString() )

    }


    def 'should increment workflow stats' () {
        given:
        def workflow = new WorkflowMetadata(
                workDir: Paths.get('/work/dir'),
                stats: new WorkflowStats(),
                nextflow: [version: '0.27.9', build: '3232', timestamp: '2017-12-12']
        )

        def file = TestHelper.createInMemTempFile('report.html')
        def observer = Spy(ReportObserver, constructorArgs: [file])
        observer.getWorkflowMetadata() >> workflow

        def TASKID1 = TaskId.of(10)
        def TASKID2 = TaskId.of(20)
        def TASKID3 = TaskId.of(30)

        def RECORD1 = new TraceRecord([task_id: TASKID1])
        def RECORD2 = new TraceRecord([task_id: TASKID2])
        def RECORD3 = new TraceRecord([task_id: TASKID3])

        def HANDLER1 = Mock(TaskHandler)
        def HANDLER2 = Mock(TaskHandler)
        def HANDLER3 = Mock(TaskHandler)

        when:
        observer.onProcessStart(HANDLER1, RECORD1)
        then:
        observer.records[TASKID1] == RECORD1

        when:
        observer.onProcessComplete(HANDLER1, RECORD1)
        then:
        1 * observer.aggregate(RECORD1) >> null

        when:
        observer.onProcessCached(HANDLER2, RECORD2)
        then:
        1 * observer.aggregate(RECORD2) >> null
        observer.records[TASKID2] == RECORD2

        when:
        observer.onProcessStart(HANDLER3, RECORD3)
        then:
        observer.records[TASKID3] == RECORD3
    }

    def 'should render not tasks payload' () {

        given:
        def observer = Spy(ReportObserver)
        def BIG = Mock(Map)
        BIG.size() >> ReportObserver.DEF_MAX_TASKS+1

        when:
        def result = observer.renderTasksJson()
        then:
        1 * observer.getRecords() >> BIG
        result == 'null'
    }

    def 'should render tasks payload' () {
        given:
        def observer = Spy(ReportObserver)

        def TASKID1 = TaskId.of(10)
        def TASKID2 = TaskId.of(20)
        def TASKID3 = TaskId.of(30)

        def RECORD1 = new TraceRecord([task_id: TASKID1, process: 'fastqc',  name:'fastqc-1', realtime: 100, '%cpu': 200, read_bytes: 300, write_bytes: 400])
        def RECORD2 = new TraceRecord([task_id: TASKID2, process: 'fastqc',  name:'fastqc-2', realtime: 500, '%cpu': 600, read_bytes: 700, write_bytes: 800])
        def RECORD3 = new TraceRecord([task_id: TASKID3, process: 'multiqc', name:'multiqc-1', realtime: 900, '%cpu': 100, read_bytes: 200, write_bytes: 300])

        def records = [
                (TASKID1): RECORD1,
                (TASKID2): RECORD2,
                (TASKID3): RECORD3
        ]


        when:
        def json = observer.renderTasksJson()
        def trace = new JsonSlurper().parseText(json) as List
        then:
        1 * observer.getRecords() >> records
        trace instanceof List
        trace.size() == 3
        trace[0].'task_id' == '10'
        trace[0].'process' == 'fastqc'
        trace[0].'name' == 'fastqc-1'
        trace[0].'realtime' == '100'
        trace[0].'%cpu' == '200'
        trace[0].'read_bytes' == '300'
        trace[0].'write_bytes' == '400'

        trace[1].'task_id' == '20'
        trace[1].'process' == 'fastqc'
        trace[1].'name' == 'fastqc-2'
        trace[1].'realtime' == '500'
        trace[1].'%cpu' == '600'
        trace[1].'read_bytes' == '700'
        trace[1].'write_bytes' == '800'

        trace[2].'task_id' == '30'
        trace[2].'process' == 'multiqc'
        trace[2].'name' == 'multiqc-1'
        trace[2].'realtime' == '900'
        trace[2].'%cpu' == '100'
        trace[2].'read_bytes' == '200'
        trace[2].'write_bytes' == '300'
    }

    def 'should render summary json' () {
        given:
        def observer = [:] as ReportObserver
        observer.executor = Executors.newCachedThreadPool()
        def r1 = new TraceRecord([process: 'bwa-mem', name: 'bwa-mem-1',          '%cpu':1_000, vmem:2_000,                   realtime:3_000,  time: 5_000,  read_bytes:4_000, write_bytes:5_000 ])
        def r2 = new TraceRecord([process: 'bwa-mem', name: 'bwa-mem-2',          '%cpu':6_000, vmem:7_000,   memory: 10_000, realtime:8_000,  time: 10_000, read_bytes:9_000, write_bytes:10_000 ])
        def r3 = new TraceRecord([process: 'bwa-mem', name: 'bwa-mem-3', cpus: 2, '%cpu':10_000, vmem:12_000, memory: 10_000, realtime:13_000, time: 10_000, read_bytes:14_000, write_bytes:15_000 ])
        def r4 = new TraceRecord([process: 'multiqc', name: 'multiqc-1',          '%cpu':16_000, vmem:17_000, memory: 20_000, realtime:18_000, time: 20_000, read_bytes:19_000, write_bytes:20_000 ])
        def r5 = new TraceRecord([process: 'multiqc', name: 'multiqc-2', cpus: 2, '%cpu':21_000, vmem:22_000, memory: 20_000, realtime:23_000, time: 20_000, read_bytes:24_000, write_bytes:25_000 ])

        observer.aggregate(r1)
        observer.aggregate(r2)
        observer.aggregate(r3)
        observer.aggregate(r4)
        observer.aggregate(r5)

        when:
        def json = observer.renderSummaryJson()
        def result = new JsonSlurper().parseText(json)
        then:
        result.size() == 2
        result.'bwa-mem' instanceof Map
        result.'bwa-mem'.cpu instanceof Map
        result.'bwa-mem'.mem instanceof Map
        result.'bwa-mem'.time instanceof Map
        result.'bwa-mem'.reads instanceof Map
        result.'bwa-mem'.writes instanceof Map

        result.'multiqc' instanceof Map
        result.'multiqc'.cpu instanceof Map
        result.'multiqc'.mem instanceof Map
        result.'multiqc'.time instanceof Map
        result.'multiqc'.reads instanceof Map
        result.'multiqc'.writes instanceof Map

        result.'bwa-mem'.time.min == 3000
        result.'bwa-mem'.time.max == 13000
        result.'bwa-mem'.time.q1 == 5500
        result.'bwa-mem'.time.q2 == 8000
        result.'bwa-mem'.time.q3 == 10500
        result.'bwa-mem'.time.mean == 8000

        result.'bwa-mem'.timeUsage.min == 60
        result.'bwa-mem'.timeUsage.max == 130
        result.'bwa-mem'.timeUsage.q1 == 70
        result.'bwa-mem'.timeUsage.q2 == 80
        result.'bwa-mem'.timeUsage.q3 == 105
        result.'bwa-mem'.timeUsage.mean == 90


        result.'bwa-mem'.cpu."min" == 1000
        result.'bwa-mem'.cpu."minLabel" == 'bwa-mem-1'
        result.'bwa-mem'.cpu."max" == 10000
        result.'bwa-mem'.cpu."maxLabel" == "bwa-mem-3"
        result.'bwa-mem'.cpu."q1" == 3500
        result.'bwa-mem'.cpu."q2" == 6000
        result.'bwa-mem'.cpu."q3" == 8000
        result.'bwa-mem'.cpu."mean" == 5666.67

        result.'bwa-mem'.cpuUsage."min" == 1000
        result.'bwa-mem'.cpuUsage."minLabel" == 'bwa-mem-1'
        result.'bwa-mem'.cpuUsage."max" == 6000
        result.'bwa-mem'.cpuUsage."maxLabel" == "bwa-mem-2"
        result.'bwa-mem'.cpuUsage."q1" == 3000
        result.'bwa-mem'.cpuUsage."q2" == 5000
        result.'bwa-mem'.cpuUsage."q3" == 5500
        result.'bwa-mem'.cpuUsage."mean" == 4000

        result.'multiqc'.mem."min" == 17000.0
        result.'multiqc'.mem."minLabel" == "multiqc-1"
        result.'multiqc'.mem."max" == 22000.0
        result.'multiqc'.mem."maxLabel" == "multiqc-2"
        result.'multiqc'.mem."q1" == 18250.0
        result.'multiqc'.mem."q2" == 19500.0
        result.'multiqc'.mem."q3" == 20750.0
        result.'multiqc'.mem."mean" == 19500.0

        result.'multiqc'.memUsage."min" == 85
        result.'multiqc'.memUsage."minLabel" == "multiqc-1"
        result.'multiqc'.memUsage."max" == 110
        result.'multiqc'.memUsage."maxLabel" == "multiqc-2"
        result.'multiqc'.memUsage."q1" == 91.25
        result.'multiqc'.memUsage."q2" == 97.50
        result.'multiqc'.memUsage."q3" == 103.75
        result.'multiqc'.memUsage."mean" == 97.50

        result.'multiqc'.time.min == 18000
        result.'multiqc'.time.max == 23000
        result.'multiqc'.time.q1 == 19250
        result.'multiqc'.time.q2 == 20500
        result.'multiqc'.time.q3 == 21750
        result.'multiqc'.time.mean == 20500

        result.'multiqc'.timeUsage.min == 90
        result.'multiqc'.timeUsage.max == 115
        result.'multiqc'.timeUsage.q1 == 96.25
        result.'multiqc'.timeUsage.q2 == 102.50
        result.'multiqc'.timeUsage.q3 == 108.75
        result.'multiqc'.timeUsage.mean == 102.50
        
        cleanup:
        observer?.executor?.shutdown()
    }
}
