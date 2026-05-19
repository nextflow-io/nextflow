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

import java.util.concurrent.CountDownLatch
import java.util.concurrent.Executors
import java.util.concurrent.TimeUnit
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

    def 'should submit gate.prepare to executor on schedule'() {
        given:
        def latch = new CountDownLatch(1)
        def gate = Mock(TaskReadinessGate) {
            1 * prepare(_) >> { latch.countDown() }
        }
        def handler = mockHandler()
        def monitor = newMonitor()
        monitor.readinessGates = [gate]
        monitor.gateExecutor = Executors.newVirtualThreadPerTaskExecutor()

        when:
        monitor.schedule(handler)

        then:
        latch.await(5, TimeUnit.SECONDS)
        monitor.gateStates.containsKey(handler)
        monitor.gateStates.get(handler).futures.size() == 1
    }

    def 'should not touch gate state when no gates are registered'() {
        given:
        def handler = mockHandler()
        def monitor = newMonitor()

        when:
        monitor.schedule(handler)

        then:
        monitor.gateStates.isEmpty()
    }

    private TaskPollingMonitor newMonitor() {
        new TaskPollingMonitor(name: 'local', session: Mock(Session), config: new ExecutorConfig([:]),
                               pollInterval: '1s', capacity: 10)
    }

    private TaskHandler mockHandler() {
        Mock(TaskHandler) { getTask() >> Mock(TaskRun) }
    }
}
