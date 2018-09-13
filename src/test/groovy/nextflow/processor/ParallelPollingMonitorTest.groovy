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

import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicInteger

import nextflow.Session
import nextflow.util.Duration
import nextflow.util.ThrottlingExecutor
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ParallelPollingMonitorTest extends Specification {

    def 'should retry task' () {

        given:
        def success = new AtomicInteger()
        def failure = new AtomicInteger()
        def count = new AtomicInteger()
        def change = new AtomicInteger()
        def retry = new AtomicInteger()

        def session = Mock(Session)
        def handler = Mock(TaskHandler)

        def opts = new ThrottlingExecutor.Options()
                        .retryOn(IllegalArgumentException)
                        .withRateLimit('10/sec')
                        .withErrorBurstDelay(Duration.of('5sec'))
                        .withAutoThrottle()
                        .onSuccess { success.incrementAndGet() }
                        .onRetry { retry.incrementAndGet() }
                        .onFailure { failure.incrementAndGet() }
                        .onRateLimitChange { change.incrementAndGet() }

        def exec = ThrottlingExecutor.create(opts)

        def reaper = Mock(ThrottlingExecutor)

        def mon = new ParallelPollingMonitor(exec, reaper, [session:session, name:'foo', pollInterval:'1sec'])
        when:
        mon.submitter = exec
        mon.submit(handler)
        sleep 3_000
        exec.shutdown()
        exec.awaitTermination(1, TimeUnit.MINUTES)

        then:
        handler.submit() >> { println "c=$count"; if(count.getAndIncrement()<2) throw new IllegalArgumentException("Ooops!") }

        success.get() == 1
        retry.get() == 2
        failure.get() == 0
        change.get() == 1

    }

}
