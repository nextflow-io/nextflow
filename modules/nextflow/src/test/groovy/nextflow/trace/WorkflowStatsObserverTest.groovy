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

import java.util.concurrent.TimeUnit

import spock.lang.Specification
import spock.lang.Timeout

@Timeout(value = 10, unit = TimeUnit.SECONDS)
class WorkflowStatsObserverTest extends Specification {

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
