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
import java.util.concurrent.TimeUnit

import nextflow.Session
import nextflow.executor.ExecutorConfig
import spock.lang.Specification

/**
 * Monitor-level integration tests for the {@link TaskGateManager} hook.
 * Behaviour of the manager itself is covered in {@link TaskGateManagerTest};
 * these tests pin only the three monitor delegation points (schedule, canSubmit, evict).
 */
class TaskPollingMonitorReadinessGateTest extends Specification {

    def 'should hold task in pending while gate runs and admit it once gate completes'() {
        given:
        def block = new CountDownLatch(1)
        def gate = Mock(TaskReadinessGate) { prepare(_) >> { block.await(5, TimeUnit.SECONDS) } }
        def handler = mockReadyHandler()
        def monitor = newMonitorWithGates([gate])

        when:
        monitor.schedule(handler)
        def busy = monitor.canSubmit(handler)
        block.countDown()
        Thread.sleep(200)
        def ready = monitor.canSubmit(handler)

        then:
        !busy
        ready
    }

    def 'should cancel in-flight gate futures on evict'() {
        given:
        def started = new CountDownLatch(1)
        def interrupted = new CountDownLatch(1)
        def gate = new TaskReadinessGate() {
            void prepare(TaskHandler h) throws InterruptedException {
                started.countDown()
                try { new CountDownLatch(1).await() }
                catch( InterruptedException e ) { interrupted.countDown(); throw e }
            }
        }
        def handler = mockReadyHandler()
        def monitor = newMonitorWithGates([gate])

        when:
        monitor.schedule(handler)
        started.await(5, TimeUnit.SECONDS)
        monitor.evict(handler)

        then:
        interrupted.await(5, TimeUnit.SECONDS)
    }

    def 'should be a zero-cost path when no gate is registered'() {
        given:
        def handler = mockReadyHandler()
        def monitor = newMonitorWithGates([])

        when:
        monitor.schedule(handler)

        then:
        monitor.canSubmit(handler)            // empty manager admits immediately
        monitor.gateManager.trackedHandlers().isEmpty()
    }

    private TaskPollingMonitor newMonitorWithGates(List<TaskReadinessGate> gates) {
        def monitor = new TaskPollingMonitor(name: 'local', session: Mock(Session), config: new ExecutorConfig([:]),
                                             pollInterval: '1s', capacity: 10)
        monitor.gateManager = new TaskGateManager(null, gates)
        return monitor
    }

    private TaskHandler mockReadyHandler() {
        Mock(TaskHandler) {
            isReady() >> true
            canForkProcess() >> true
            getTask() >> Mock(TaskRun) { getName() >> 'task-x' }
        }
    }
}
