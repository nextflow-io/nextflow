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

package nextflow.executor

import java.util.concurrent.CancellationException
import java.util.concurrent.CountDownLatch
import java.util.concurrent.CyclicBarrier
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicInteger

import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import spock.lang.Specification
import spock.lang.Timeout
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Timeout(30)
class ArrayLauncherPoolTest extends Specification {

    def 'should prepare launchers in parallel' () {
        given:
        def pool = new ArrayLauncherPool('test', 4)
        def barrier = new CyclicBarrier(4)
        def handlers = (1..4).collect { Mock(TaskHandler) { prepareLauncher() >> { barrier.await(10, TimeUnit.SECONDS) } } }

        when:
        pool.prepare(handlers)
        then:
        noExceptionThrown()

        cleanup:
        pool.shutdown()
    }

    def 'should bound the number of concurrent preparations' () {
        given:
        def pool = new ArrayLauncherPool('test', 2)
        def running = new AtomicInteger()
        def peak = new AtomicInteger()
        def handlers = (1..10).collect {
            Mock(TaskHandler) {
                prepareLauncher() >> {
                    peak.accumulateAndGet(running.incrementAndGet(), Math::max)
                    sleep 50
                    running.decrementAndGet()
                }
            }
        }

        when:
        pool.prepare(handlers)
        then:
        peak.get() == 2

        cleanup:
        pool.shutdown()
    }

    def 'should fail as soon as a child fails' () {
        given:
        def pool = new ArrayLauncherPool('test', 2)
        def blocked = new CountDownLatch(1)
        def slow = Mock(TaskHandler) { prepareLauncher() >> { blocked.await() } }
        def fail = Mock(TaskHandler) {
            prepareLauncher() >> { throw new IOException('boom') }
            getTask() >> Mock(TaskRun)
        }

        when:
        pool.prepare([slow, fail])
        then:
        def e = thrown(IOException)
        e.message == 'boom'

        cleanup:
        pool.shutdown()
    }

    def 'should not block the caller when the pool is shut down' () {
        given:
        def pool = new ArrayLauncherPool('test', 1)
        def started = new CountDownLatch(1)
        def slow = Mock(TaskHandler) { prepareLauncher() >> { started.countDown(); sleep 60_000 } }
        def queued = Mock(TaskHandler)
        Throwable error = null
        def caller = Thread.start {
            try { pool.prepare([slow, queued]) }
            catch( Throwable t ) { error = t }
        }

        when:
        started.await()
        pool.shutdown()
        caller.join(10_000)
        then:
        !caller.isAlive()
        error instanceof CancellationException || error instanceof InterruptedException
    }

}
