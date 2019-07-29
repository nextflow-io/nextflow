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
import java.nio.file.Paths
import java.time.OffsetDateTime

import groovy.json.JsonSlurper
import nextflow.NextflowMeta
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
                start: OffsetDateTime.now(),
                complete: OffsetDateTime.now(),
                workDir: workDir,
                stats: new WorkflowStats(),
                nextflow: new NextflowMeta('0.27.9', 3232, '2017-12-12')
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
                nextflow: new NextflowMeta('0.27.9', 3232, '2017-12-12')
        )

        def aggregator = Mock(ResourcesAggregator)
        def file = TestHelper.createInMemTempFile('report.html')
        ReportObserver observer = Spy(ReportObserver, constructorArgs: [file])
        observer.getWorkflowMetadata() >> workflow
        observer.@aggregator = aggregator

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


}
