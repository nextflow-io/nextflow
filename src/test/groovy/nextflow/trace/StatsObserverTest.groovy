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
        observer.onProcessComplete(HANDLER1)
        then:
        1 * HANDLER1.getTraceRecord() >> RECORD1
        1 * stats.updateTasksCompleted(RECORD1) >> null

        when:
        observer.onProcessCached(HANDLER2)
        then:
        1 * HANDLER2.getTraceRecord() >> RECORD2
        1 * stats.updateTasksCached(RECORD2) >> null
    }


}
