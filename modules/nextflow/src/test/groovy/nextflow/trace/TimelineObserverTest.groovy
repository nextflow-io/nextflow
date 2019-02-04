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
import java.nio.file.Path
import java.nio.file.Paths

import nextflow.processor.TaskHandler
import nextflow.processor.TaskId
import spock.lang.Specification
import test.TestHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TimelineObserverTest extends Specification {

    def 'should read html template' () {

        given:
        def observer = [:] as TimelineObserver
        when:
        def tpl = observer.readTemplate()
        then:
        tpl.startsWith '<!doctype html>'
        tpl.contains 'REPLACE_WITH_TIMELINE_DATA'

    }

    def 'should return timeline in json format' () {

        given:
        def now = 1429821425141
        def r1 = new TraceRecord()
        r1.task_id = '1'
        r1.name = 'foo'
        r1.process = 'alpha'
        r1.peak_rss = 50_000_000

        def r2 = new TraceRecord()
        r2.task_id = '2'
        r2.name = 'bar'
        r2.submit = now
        r2.start = now + 100
        r2.complete = now + 500
        r2.realtime = 400
        r2.duration = 500
        r2.process = 'alpha'
        r2.peak_rss = 60_000_000

        def r3 = new TraceRecord()
        r3.task_id = '3'
        r3.name = 'baz'
        r3.submit = now
        r3.start = now + 200
        r3.complete = now + 700
        r3.realtime = 500
        r3.duration = 700
        r3.process = 'beta'
        r3.peak_rss = 70_000_000

        def observer = [:] as TimelineObserver
        observer.beginMillis = 1000
        observer.startMillis = 1000
        observer.endMillis = 3500
        observer.records['1'] = r1
        observer.records['2'] = r2
        observer.records['3'] = r3

        expect:
        observer.renderData().toString() == /
            var elapsed="2.5s"
            var beginningMillis=1000;
            var endingMillis=3500;
            var data=[
            {"label": "foo", "times": []},
            {"label": "bar", "times": [{"starting_time": 1429821425141, "ending_time": 1429821425241, "color":c1(0)}, {"starting_time": 1429821425241, "ending_time": 1429821425641, "color":c2(0), "label": "500ms \\/ 57.2 MB"}]},
            {"label": "baz", "times": [{"starting_time": 1429821425141, "ending_time": 1429821425341, "color":c1(1)}, {"starting_time": 1429821425341, "ending_time": 1429821425841, "color":c2(1), "label": "700ms \\/ 66.8 MB"}]}
            ]
            /
            .stripIndent().leftTrim()
    }

    def 'should add records' () {

        given:
        def now = 1429821425141
        def r1 = new TraceRecord()
        r1.task_id = TaskId.of(1)
        r1.name = 'foo'
        r1.process = 'alpha'

        def r2 = new TraceRecord()
        r2.task_id = TaskId.of(2)
        r2.name = 'bar'
        r2.submit = now
        r2.start = now + 100
        r2.complete = now + 500
        r2.realtime = 400
        r2.duration = 500
        r2.process = 'alpha'

        def r3 = new TraceRecord()
        r3.task_id = TaskId.of(3)
        r3.name = 'baz'
        r3.submit = now
        r3.start = now + 200
        r3.complete = now + 700
        r3.realtime = 500
        r3.duration = 700
        r3.process = 'beta'

        def h1 = Mock(TaskHandler)
        h1.getTask() >> [id: TaskId.of(1)]
        h1.getTraceRecord() >> r1

        def h2 = Mock(TaskHandler)
        h2.getTask() >> [id: TaskId.of(2)]
        h2.getTraceRecord() >> r2

        def h3 = Mock(TaskHandler)
        h3.getTask() >> [id: TaskId.of(3)]
        h3.getTraceRecord() >> r3

        when:
        def observer = new TimelineObserver(Mock(Path))
        observer.onProcessComplete(h1, h1.getTraceRecord())
        observer.onProcessComplete(h2, h2.getTraceRecord())
        observer.onProcessComplete(h3, h3.getTraceRecord())
        then:
        observer.records[TaskId.of(1)] == r1
        observer.records[TaskId.of(2)] == r2
        observer.records[TaskId.of(3)] == r3

    }

    def 'should create html file' () {

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
        r2.peak_rss = 60_000_000

        def r3 = new TraceRecord()
        r3.task_id = '3'
        r3.name = 'baz'
        r3.submit = now
        r3.start = now + 200
        r3.complete = now + 700
        r3.realtime = 500
        r3.duration = 700
        r3.process = 'beta'
        r3.peak_rss = 70_000_000

        def file = TestHelper.createInMemTempFile('report.html')
        def observer = new TimelineObserver(file)
        observer.beginMillis = 1000
        observer.startMillis = 1000
        observer.endMillis = 3500
        observer.records['1'] = r1
        observer.records['2'] = r2
        observer.records['3'] = r3

        when:
        observer.renderHtml()
        then:
        Files.exists(file)
        file.text == Paths.get('src/test/groovy/nextflow/trace/timeline-expected.html').text

    }

}
