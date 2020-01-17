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
import nextflow.processor.TaskProcessor
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class StatsObserverTest extends Specification {

    def 'should create process' () {
        given:
        def observer = new WorkflowStatsObserver(Mock(Session))
        and:
        def process = Mock(TaskProcessor) { getId() >> 1; getName() >> 'foo' }

        when:
        observer.onProcessCreate(process)
        then:
        observer.stats.processes.size() == 1
        and:
        with(observer.stats.processes[0]) {
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
        def observer = new WorkflowStatsObserver(Mock(Session))
        def ts = System.currentTimeMillis()
        and:
        def process = Mock(TaskProcessor) { getId() >> 1; getName() >> 'foo' }
        and:
        observer.onProcessCreate(process)

        when:
        observer.onProcessTerminate(process)
        then:
        with(observer.stats.processes[0]) {
            name == 'foo'
            terminated
        }
        and:
        observer.getChangeTimestamp() >= ts
        observer.getChangeTimestamp() <= System.currentTimeMillis()
    }


}
