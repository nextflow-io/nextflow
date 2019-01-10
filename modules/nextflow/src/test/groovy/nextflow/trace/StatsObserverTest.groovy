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
import nextflow.processor.TaskHandler
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class StatsObserverTest extends Specification {

    def 'should update tasks completed' () {
        given:
        def stats = Mock(WorkflowStats)
        def observer = new StatsObserver(stats: stats)

        def RECORD1 = Mock(TraceRecord)
        def RECORD2 = Mock(TraceRecord)

        def HANDLER1 = Mock(TaskHandler)
        def HANDLER2 = Mock(TaskHandler)

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
