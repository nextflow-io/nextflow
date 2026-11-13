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

package nextflow.trace

import java.util.concurrent.CountDownLatch
import java.util.concurrent.TimeUnit

import nextflow.Session
import nextflow.processor.TaskProcessor
import spock.lang.Specification
import spock.lang.Timeout

@Timeout(value = 10, unit = TimeUnit.SECONDS)
class WorkflowStatsObserverTest extends Specification {

    private TaskProcessor process(int id) {
        Mock(TaskProcessor) { getId() >> id; getName() >> "proc_$id" }
    }

    private TaskProcessor failingProcess() {
        Mock(TaskProcessor) { getId() >> 0; getName() >> { throw new IllegalStateException('boom') } }
    }

    def 'should apply the updates sent from several threads' () {
        given:
        def observer = new WorkflowStatsObserver(Mock(Session))
        def processes = (1..500).collect { process(it) }

        when:
        def threads = (0..<4).collect { t ->
            Thread.start {
                for( int i = t; i < processes.size(); i += 4 ) {
                    observer.onProcessCreate(processes[i])
                    observer.onProcessTerminate(processes[i])
                }
            }
        }
        threads*.join()
        def stats = observer.getStats()

        then:
        stats.processes.size() == 500
        stats.processes.every { it.terminated }
    }

    def 'should return copies that are not changed by later updates' () {
        given:
        def observer = new WorkflowStatsObserver(Mock(Session))
        def foo = process(1)
        observer.onProcessCreate(foo)

        when:
        def stats = observer.getStats()
        def quick = observer.getQuickStats()
        observer.onProcessTerminate(foo)
        observer.onProcessCreate(process(2))
        def latest = observer.getStats()

        then:
        stats.processes.size() == 1
        !stats.processes[0].terminated
        quick.processes.size() == 1
        !quick.processes[0].terminated
        and:
        latest.processes.size() == 2
        latest.processes.find { it.name == 'proc_1' }.terminated
    }

    def 'should abort the session when an update fails' () {
        given:
        def aborted = new CountDownLatch(1)
        def session = Mock(Session)
        def observer = new WorkflowStatsObserver(session)

        when:
        observer.onProcessCreate(failingProcess())
        aborted.await(5, TimeUnit.SECONDS)

        then:
        1 * session.abort({ it instanceof IllegalStateException }) >> { aborted.countDown() }
    }

    def 'should return the stats when the abort waits for them' () {
        given:
        WorkflowStatsObserver observer
        WorkflowStats direct
        WorkflowStats joined
        def done = new CountDownLatch(1)
        def session = Mock(Session) {
            abort(_) >> {
                // read the stats on the aborting thread, and on another thread joined by it
                direct = observer.getStats()
                def reader = Thread.start { joined = observer.getStats() }
                reader.join()
                done.countDown()
            }
        }
        observer = new WorkflowStatsObserver(session)

        when:
        observer.onProcessCreate(process(1))
        observer.onProcessCreate(failingProcess())

        then:
        done.await(5, TimeUnit.SECONDS)
        direct.processes.size() == 1
        joined.processes.size() == 1
    }

    def 'should return the quick stats and keep the interrupt when interrupted' () {
        given:
        def observer = new WorkflowStatsObserver(null)
        def stats
        def interrupted

        when:
        Thread.currentThread().interrupt()
        try {
            stats = observer.getStats()
            interrupted = Thread.currentThread().isInterrupted()
        }
        finally {
            Thread.interrupted()
        }

        then:
        stats != null
        interrupted
        observer.getStats() != null
    }
}
