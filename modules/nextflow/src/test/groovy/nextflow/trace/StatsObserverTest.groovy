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

import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class StatsObserverTest extends Specification {

    def 'should create process' () {
        given:
        def observer = new StatsObserver(Mock(Session))
        and:
        def process = Mock(TaskProcessor) { getId() >> 1; getName() >> 'foo' }

        when:
        observer.onProcessCreate(process)
        then:
        observer.getProgressLength() == 1
        with(observer.getProgress()[0]) {
            name == 'foo'
            pending == 0
            running == 0
            submitted == 0
            failed == 0
            !terminated
        }
    }

    def 'should terminate process' () {
        given:
        def observer = new StatsObserver(Mock(Session))
        def ts = System.currentTimeMillis()
        and:
        def process = Mock(TaskProcessor) { getId() >> 1; getName() >> 'foo' }
        and:
        observer.onProcessCreate(process)

        when:
        observer.onProcessTerminate(process)
        then:
        observer.getProgressLength() == 1
        with(observer.getProgress()[0]) {
            name == 'foo'
            terminated
        }
        and:
        observer.getChangeTimestamp() >= ts
        observer.getChangeTimestamp() <= System.currentTimeMillis()
    }

    def 'should mark pending' () {
        given:
        def observer = new StatsObserver(Mock(Session))
        and:
        def task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) {
                getId() >> 1;
                getName() >> 'foo' }
        }
        and:
        def handler = Mock(TaskHandler) { getTask() >> task }
        and:
        def record = Mock(ProgressRecord)
        observer.records = [ record ]

        when:
        observer.onProcessPending(handler, null)
        then:
        1 * record.markPending() >> null
    }

    def 'should mark submit' () {
        given:
        def observer = new StatsObserver(Mock(Session))
        and:
        def task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) {
                getId() >> 1;
                getName() >> 'foo' }
        }
        and:
        def handler = Mock(TaskHandler) { getTask() >> task }
        and:
        def record = Mock(ProgressRecord)
        observer.records = [ record ]

        when:
        observer.onProcessSubmit(handler, null)
        then:
        1 * record.markSubmitted(task) >> null
    }

    def 'should mark running' () {
        given:
        def observer = new StatsObserver(Mock(Session))
        and:
        def task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) {
                getId() >> 1;
                getName() >> 'foo' }
        }
        and:
        def handler = Mock(TaskHandler) { getTask() >> task }
        and:
        def record = Mock(ProgressRecord)
        observer.records = [ record ]

        when:
        observer.onProcessStart(handler, null)
        then:
        1 * record.markRunning(task) >> null
    }

    def 'should mark complete' () {
        given:
        def observer = new StatsObserver(Mock(Session))
        and:
        def task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) {
                getId() >> 1;
                getName() >> 'foo' }
        }
        and:
        def handler = Mock(TaskHandler) { getTask() >> task }
        and:
        def record = Mock(ProgressRecord)
        observer.records = [ record ]

        when:
        observer.onProcessComplete(handler, null)
        then:
        1 * record.markComplete(task) >> null
    }


    def 'should update tasks stats' () {
        given:
        def stats = Mock(WorkflowStats)
        def observer = new StatsObserver(Mock(Session))
        observer.stats = stats

        def RECORD1 = Mock(TraceRecord)
        def RECORD2 = Mock(TraceRecord)

        def HANDLER1 = Mock(TaskHandler) {
            getTask() >> Mock(TaskRun) {
                getProcessor() >> Mock(TaskProcessor) { getId() >> 1 }
            }
        }
        def HANDLER2 = Mock(TaskHandler) {
            getTask() >> Mock(TaskRun) {
                getProcessor() >> Mock(TaskProcessor) { getId() >> 1 }
            }
        }

        when:
        observer.onProcessComplete(HANDLER1, RECORD1)
        then:
        1 * stats.updateTasksCompleted(RECORD1) >> null

        when:
        observer.onProcessCached(HANDLER2, RECORD2)
        then:
        1 * stats.updateTasksCached(RECORD2) >> null
    }


}
