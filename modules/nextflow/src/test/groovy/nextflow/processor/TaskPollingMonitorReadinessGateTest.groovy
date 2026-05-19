/*
 * Copyright 2013-2026, Seqera Labs
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

import java.util.concurrent.Future
import nextflow.Session
import nextflow.executor.ExecutorConfig
import spock.lang.Specification

class TaskPollingMonitorReadinessGateTest extends Specification {

    def 'should leave gateExecutor null when no gates are registered'() {
        given:
        def session = Mock(Session)
        def config = new ExecutorConfig([:])
        def monitor = new TaskPollingMonitor(name: 'local', session: session, config: config, pollInterval: '1s', capacity: 10)

        when:
        monitor.readinessGates = []   // simulate empty extension list

        then:
        monitor.gateExecutor == null
        monitor.gateStates.isEmpty()
    }

    def 'should expose GateState as a package-scoped type'() {
        when:
        def state = new TaskPollingMonitor.GateState([])
        then:
        state.futures == []
        state.scheduledAt > 0
    }
}
