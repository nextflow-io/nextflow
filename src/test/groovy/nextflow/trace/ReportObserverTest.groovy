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

package nextflow.trace
import java.nio.file.Files
import java.nio.file.Paths

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
        def json = new JsonSlurper().parseText(str)

        then:
        json.trace[0].task_id == '1'
        json.trace[0].name == 'foo'
        json.trace[0].process == 'alpha'

        json.trace[1].task_id == '2'
        json.trace[1].name == 'bar'
        json.trace[1].submit == '1429821425141'
        json.trace[1].start == '1429821425241'
        json.trace[1].complete == '1429821425641'
        json.trace[1].realtime == '400'
        json.trace[1].duration == '500'
        json.trace[1].process == 'alpha'

        json.trace[2].task_id == '3'
        json.trace[2].name == 'baz'
        json.trace[2].submit == '1429821425141'
        json.trace[2].start == '1429821425341'
        json.trace[2].complete == '1429821425841'
        json.trace[2].realtime == '500'
        json.trace[2].duration == '700'
        json.trace[2].process == 'beta'

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
        observer.onProcessStart(HANDLER1)
        then:
        1 * HANDLER1.getTraceRecord() >> RECORD1
        observer.records[TASKID1] == RECORD1

        when:
        observer.onProcessComplete(HANDLER1)
        then:
        1 * HANDLER1.getTraceRecord() >> RECORD1

        when:
        observer.onProcessCached(HANDLER2)
        then:
        1 * HANDLER2.getTraceRecord() >> RECORD2
        observer.records[TASKID2] == RECORD2

        when:
        observer.onProcessStart(HANDLER3)
        then:
        1 * HANDLER3.getTraceRecord() >> RECORD3
        observer.records[TASKID3] == RECORD3
    }
}
