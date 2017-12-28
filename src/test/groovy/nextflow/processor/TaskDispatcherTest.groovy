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

package nextflow.processor

import nextflow.executor.SgeExecutor
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskDispatcherTest extends Specification {

    def testGetMonitor()  {

        setup:
        def dispatcher = [:] as TaskDispatcher
        def monitor1 = Mock(TaskPollingMonitor)
        def monitor2 = Mock(TaskPollingMonitor)
        def executor = [:] as SgeExecutor

        expect:
        monitor1 == dispatcher.getOrCreateMonitor(SgeExecutor) { monitor1 }
        monitor2 != dispatcher.getOrCreateMonitor(SgeExecutor) { monitor1 }
        monitor1 == dispatcher.getMonitor(SgeExecutor)
        monitor1 == dispatcher.getMonitor(executor)

    }

}
