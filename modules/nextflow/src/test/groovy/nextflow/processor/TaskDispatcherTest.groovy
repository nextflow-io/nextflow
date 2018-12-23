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
