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
import nextflow.exception.ProcessException
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

    def 'should return false from canSubmit while gate prepare runs'() {
        given:
        def block = new CountDownLatch(1)
        def gate = Mock(TaskReadinessGate) {
            prepare(_) >> { block.await(5, TimeUnit.SECONDS) }
        }
        def handler = mockReadyHandler()
        def monitor = newMonitor()
        monitor.readinessGates = [gate]
        monitor.gateExecutor = Executors.newVirtualThreadPerTaskExecutor()

        when:
        monitor.schedule(handler)
        def busy = monitor.canSubmit(handler)
        block.countDown()
        Thread.sleep(200)
        def ready = monitor.canSubmit(handler)

        then:
        !busy
        ready
        !monitor.gateStates.containsKey(handler)   // state removed on admission
    }

    def 'should rethrow ProcessException from gate as-is'() {
        given:
        // ProcessException (and subclasses) propagate unchanged so that downstream
        // code in TaskProcessor.resumeOrDie sees the original type and message.
        def boom = new ProcessException('boom')
        def gate = new TaskReadinessGate() {
            void prepare(TaskHandler h) throws InterruptedException {
                throw boom
            }
        }
        def handler = mockReadyHandler()
        def monitor = newMonitor()
        monitor.readinessGates = [gate]
        monitor.gateExecutor = Executors.newVirtualThreadPerTaskExecutor()

        when:
        monitor.schedule(handler)
        awaitFutureCompletion(monitor, handler)
        monitor.canSubmit(handler)

        then:
        def e = thrown(ProcessException)
        e.is(boom)                                // same instance, not wrapped
        !monitor.gateStates.containsKey(handler)  // cleaned up
    }

    def 'should preserve arbitrary RuntimeException from gate without wrapping'() {
        given:
        def gate = Mock(TaskReadinessGate) {
            prepare(_) >> { throw new IllegalStateException('underlying') }
        }
        def handler = mockReadyHandler()
        def monitor = newMonitor()
        monitor.readinessGates = [gate]
        monitor.gateExecutor = Executors.newVirtualThreadPerTaskExecutor()

        when:
        monitor.schedule(handler)
        awaitFutureCompletion(monitor, handler)
        monitor.canSubmit(handler)

        then:
        def e = thrown(IllegalStateException)
        e.message == 'underlying'
        !monitor.gateStates.containsKey(handler)
    }

    def 'should wrap checked exception causes in a ProcessException'() {
        given:
        def gate = new TaskReadinessGate() {
            void prepare(TaskHandler h) throws InterruptedException {
                throw new IOException('disk failed')
            }
        }
        def handler = mockReadyHandler()
        def monitor = newMonitor()
        monitor.readinessGates = [gate]
        monitor.gateExecutor = Executors.newVirtualThreadPerTaskExecutor()

        when:
        monitor.schedule(handler)
        awaitFutureCompletion(monitor, handler)
        monitor.canSubmit(handler)

        then:
        def e = thrown(ProcessException)
        e.message.contains('Task readiness gate failed')
        e.cause instanceof IOException
        e.cause.message == 'disk failed'
        !monitor.gateStates.containsKey(handler)
    }

    def 'should throw ProcessException when a gate future is cancelled externally'() {
        given:
        def block = new CountDownLatch(1)
        def gate = Mock(TaskReadinessGate) {
            prepare(_) >> { block.await(5, TimeUnit.SECONDS) }
        }
        def handler = mockReadyHandler()
        def monitor = newMonitor()
        monitor.readinessGates = [gate]
        monitor.gateExecutor = Executors.newVirtualThreadPerTaskExecutor()

        when:
        monitor.schedule(handler)
        monitor.gateStates.get(handler).futures.each { it.cancel(true) }
        Thread.sleep(50)
        monitor.canSubmit(handler)

        then:
        def e = thrown(ProcessException)
        e.message.contains('cancelled')
        !monitor.gateStates.containsKey(handler)

        cleanup:
        block.countDown()
    }

    // Note on ExecutionException with null cause: the allGatesReady code path
    // (lines ~319-321 in TaskPollingMonitor) wraps a null-cause ExecutionException
    // in a ProcessException("Task readiness gate failed...", e). This branch is
    // reachable but awkward to trigger from a real gate Mock — the JDK's
    // FutureTask.setException(null) throws NullPointerException, so we can't
    // surface a null-cause ExecutionException without bending the test. The
    // branch is verified by code-read.

    private void awaitFutureCompletion(TaskPollingMonitor monitor, TaskHandler handler) {
        final state = monitor.gateStates.get(handler)
        for( f in state.futures ) {
            try { f.get(5, TimeUnit.SECONDS) }
            catch( Exception ignored ) { /* swallow — we only want to wait until done */ }
        }
    }

    private TaskPollingMonitor newMonitor() {
        new TaskPollingMonitor(name: 'local', session: Mock(Session), config: new ExecutorConfig([:]),
                               pollInterval: '1s', capacity: 10)
    }

    private TaskHandler mockHandler() {
        Mock(TaskHandler) { getTask() >> Mock(TaskRun) }
    }

    private TaskHandler mockReadyHandler() {
        Mock(TaskHandler) {
            isReady() >> true
            canForkProcess() >> true
            getTask() >> Mock(TaskRun) { getName() >> 'task-x' }
        }
    }
}
