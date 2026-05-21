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
import java.util.concurrent.ExecutionException
import java.util.concurrent.TimeUnit

import nextflow.Session
import nextflow.exception.ProcessException
import nextflow.exception.ProcessRetryableException
import spock.lang.Specification

class TaskGateManagerTest extends Specification {

    def 'should be a no-op when no gates are registered'() {
        given:
        def manager = new TaskGateManager(null, [])
        def handler = mockHandler()

        when:
        manager.prepare(handler)

        then:
        manager.isReady(handler)            // no gates => always ready
        manager.trackedHandlers().isEmpty() // no state recorded
        manager.evict(handler)              // safe no-op
        noExceptionThrown()
    }

    def 'should submit gate.prepare to the executor on prepare'() {
        given:
        def latch = new CountDownLatch(1)
        def gate = Mock(TaskReadinessGate) {
            1 * prepare(_) >> { latch.countDown() }
        }
        def manager = new TaskGateManager(null, [gate])
        def handler = mockHandler()

        when:
        manager.prepare(handler)

        then:
        latch.await(5, TimeUnit.SECONDS)
        manager.futuresFor(handler).size() == 1
    }

    def 'should return false from isReady while gate prepare runs'() {
        given:
        def block = new CountDownLatch(1)
        def gate = Mock(TaskReadinessGate) {
            prepare(_) >> { block.await(5, TimeUnit.SECONDS) }
        }
        def manager = new TaskGateManager(null, [gate])
        def handler = mockHandler()

        when:
        manager.prepare(handler)
        def busy = manager.isReady(handler)
        block.countDown()
        awaitFutureCompletion(manager, handler)
        def ready = manager.isReady(handler)

        then:
        !busy
        ready
        !manager.trackedHandlers().contains(handler)   // state removed on admission
    }

    def 'should rethrow ProcessException from gate as-is'() {
        given:
        def boom = new ProcessException('boom')
        def gate = new TaskReadinessGate() {
            void prepare(TaskHandler h) throws InterruptedException {
                throw boom
            }
        }
        def manager = new TaskGateManager(null, [gate])
        def handler = mockHandler()

        when:
        manager.prepare(handler)
        awaitFutureCompletion(manager, handler)
        manager.isReady(handler)

        then:
        def e = thrown(ProcessException)
        e.is(boom)
        !manager.trackedHandlers().contains(handler)
    }

    def 'should wrap RuntimeException causes in a ProcessException with the original as cause'() {
        given:
        def gate = Mock(TaskReadinessGate) {
            prepare(_) >> { throw new IllegalStateException('underlying') }
        }
        def manager = new TaskGateManager(null, [gate])
        def handler = mockHandler()

        when:
        manager.prepare(handler)
        awaitFutureCompletion(manager, handler)
        manager.isReady(handler)

        then:
        def e = thrown(ProcessException)
        e.message.contains('Task readiness gate failed')
        e.cause instanceof IllegalStateException
        e.cause.message == 'underlying'
        !manager.trackedHandlers().contains(handler)
    }

    def 'should expose ProcessRetryableException markers via error.cause for resumeOrDie'() {
        given:
        def marker = new RetryableBoom()
        def gate = Mock(TaskReadinessGate) { prepare(_) >> { throw marker } }
        def manager = new TaskGateManager(null, [gate])
        def handler = mockHandler()

        when:
        manager.prepare(handler)
        awaitFutureCompletion(manager, handler)
        manager.isReady(handler)

        then:
        def e = thrown(ProcessException)
        e.cause.is(marker)
        e.cause instanceof ProcessRetryableException
    }

    def 'should wrap checked exception causes in a ProcessException'() {
        given:
        def gate = new TaskReadinessGate() {
            void prepare(TaskHandler h) throws InterruptedException {
                throw new IOException('disk failed')
            }
        }
        def manager = new TaskGateManager(null, [gate])
        def handler = mockHandler()

        when:
        manager.prepare(handler)
        awaitFutureCompletion(manager, handler)
        manager.isReady(handler)

        then:
        def e = thrown(ProcessException)
        e.message.contains('Task readiness gate failed')
        e.cause instanceof IOException
        e.cause.message == 'disk failed'
        !manager.trackedHandlers().contains(handler)
    }

    def 'should throw ProcessException when a gate future is cancelled externally'() {
        given:
        def block = new CountDownLatch(1)
        def gate = Mock(TaskReadinessGate) {
            prepare(_) >> { block.await(5, TimeUnit.SECONDS) }
        }
        def manager = new TaskGateManager(null, [gate])
        def handler = mockHandler()

        when:
        manager.prepare(handler)
        manager.futuresFor(handler).each { it.cancel(true) }
        awaitFutureCompletion(manager, handler)
        manager.isReady(handler)

        then:
        def e = thrown(ProcessException)
        e.message.contains('cancelled')
        !manager.trackedHandlers().contains(handler)

        cleanup:
        block.countDown()
    }

    def 'should cancel peer gate futures when one gate throws'() {
        given:
        def peerStarted = new CountDownLatch(1)
        def peerInterrupted = new CountDownLatch(1)
        def failing = new TaskReadinessGate() {
            void prepare(TaskHandler h) throws InterruptedException {
                throw new ProcessException('failing gate')
            }
        }
        def peer = new TaskReadinessGate() {
            void prepare(TaskHandler h) throws InterruptedException {
                peerStarted.countDown()
                try { new CountDownLatch(1).await() }
                catch( InterruptedException e ) { peerInterrupted.countDown(); throw e }
            }
        }
        def manager = new TaskGateManager(null, [failing, peer])
        def handler = mockHandler()

        when:
        manager.prepare(handler)
        peerStarted.await(5, TimeUnit.SECONDS)
        try { manager.futuresFor(handler)[0].get(5, TimeUnit.SECONDS) }
        catch( ExecutionException ignored ) { /* expected — failing gate threw */ }

        then:
        try { manager.isReady(handler) } catch( ProcessException expected ) { /* expected */ }
        peerInterrupted.await(5, TimeUnit.SECONDS)
        !manager.trackedHandlers().contains(handler)
    }

    def 'should admit task only after every gate has completed'() {
        given:
        def latch1 = new CountDownLatch(1)
        def latch2 = new CountDownLatch(1)
        def g1 = Mock(TaskReadinessGate) { prepare(_) >> { latch1.await(5, TimeUnit.SECONDS) } }
        def g2 = Mock(TaskReadinessGate) { prepare(_) >> { latch2.await(5, TimeUnit.SECONDS) } }
        def manager = new TaskGateManager(null, [g1, g2])
        def handler = mockHandler()

        when:
        manager.prepare(handler)
        Thread.sleep(50)
        def neither = manager.isReady(handler)
        latch1.countDown()
        Thread.sleep(50)
        def onlyFirst = manager.isReady(handler)
        latch2.countDown()
        Thread.sleep(50)
        def both = manager.isReady(handler)

        then:
        !neither
        !onlyFirst
        both
        !manager.trackedHandlers().contains(handler)
    }

    def 'should cancel gate futures on evict'() {
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
        def manager = new TaskGateManager(null, [gate])
        def handler = mockHandler()

        when:
        manager.prepare(handler)
        started.await(5, TimeUnit.SECONDS)
        manager.evict(handler)

        then:
        interrupted.await(5, TimeUnit.SECONDS)
        !manager.trackedHandlers().contains(handler)
    }

    /**
     * Marker class used to verify that ProcessRetryableException-carrying causes
     * reach resumeOrDie via the outer ProcessException's `cause` field.
     */
    static class RetryableBoom extends RuntimeException implements ProcessRetryableException {
        RetryableBoom() { super('retry me') }
    }

    private TaskHandler mockHandler() {
        Mock(TaskHandler) { getTask() >> Mock(TaskRun) { getName() >> 'task-x' } }
    }

    private void awaitFutureCompletion(TaskGateManager manager, TaskHandler handler) {
        final futures = manager.futuresFor(handler)
        if( !futures ) return
        for( f in futures ) {
            try { f.get(5, TimeUnit.SECONDS) }
            catch( Exception ignored ) { /* swallow — we only want to wait until done */ }
        }
    }
}
