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

import nextflow.exception.ProcessUnrecoverableException
import nextflow.util.MemoryUnit
import spock.lang.Specification

class ResourceTrackerTest extends Specification {

    private TaskHandler handler(int cpus, String memory) {
        Mock(TaskHandler) {
            getTask() >> new TaskRun(config: new TaskConfig(cpus: cpus, memory: memory ? MemoryUnit.of(memory) : null))
        }
    }

    def 'should acquire and release resources'() {
        given:
        def tracker = new ResourceTracker(8, MemoryUnit.of('16GB').toBytes())
        def h = handler(4, '8GB')

        when:
        tracker.acquire(h)
        then:
        tracker.availableCpus() == 4
        tracker.availableMemory() == MemoryUnit.of('8GB').toBytes()
        tracker.canAcquire(h)

        when:
        tracker.acquire(h)
        then:
        tracker.availableCpus() == 0
        tracker.availableMemory() == 0
        !tracker.canAcquire(h)
        !tracker.canAcquire(handler(1, '1GB'))

        when:
        tracker.release(h)
        then:
        tracker.availableCpus() == 4
        tracker.availableMemory() == MemoryUnit.of('8GB').toBytes()
        tracker.canAcquire(h)
    }

    def 'should limit only cpus or only memory'() {
        given:
        def cpusOnly = new ResourceTracker(4, 0)
        def memoryOnly = new ResourceTracker(0, MemoryUnit.of('16GB').toBytes())

        when:
        cpusOnly.acquire(handler(2, '8GB'))
        memoryOnly.acquire(handler(2, '8GB'))
        then:
        cpusOnly.availableCpus() == 2
        cpusOnly.availableMemory() == 0
        memoryOnly.availableCpus() == 0
        memoryOnly.availableMemory() == MemoryUnit.of('8GB').toBytes()

        expect:
        new ResourceTracker(4, 0).canAcquire(handler(2, '1TB'))
        !new ResourceTracker(4, 0).canAcquire(handler(5, '1GB'))
        new ResourceTracker(0, MemoryUnit.of('16GB').toBytes()).canAcquire(handler(100, '8GB'))
        !new ResourceTracker(0, MemoryUnit.of('16GB').toBytes()).canAcquire(handler(1, '17GB'))
    }

    def 'should fail when a task exceeds the total'() {
        given:
        def tracker = new ResourceTracker(8, MemoryUnit.of('16GB').toBytes())

        when:
        tracker.validate(handler(CPUS, MEMORY))
        then:
        def e = thrown(ProcessUnrecoverableException)
        e.message == EXPECTED

        where:
        CPUS | MEMORY  | EXPECTED
        10   | '8GB'   | 'Task requirement exceeds available CPUs -- req: 10; avail: 8'
        4    | '20GB'  | 'Task requirement exceeds available memory -- req: 20 GB; avail: 16 GB'
    }

    def 'should fail when a job array exceeds the total'() {
        given:
        def tracker = new ResourceTracker(32, 0)
        def children = (1..100).collect { handler(1, null) }
        def array = Mock(TaskHandler) {
            getTask() >> Mock(TaskArrayRun) { getChildren() >> children }
        }

        when:
        tracker.validate(array)
        then:
        def e = thrown(ProcessUnrecoverableException)
        e.message == 'Array requirement exceeds available CPUs -- req: 100 (array size: 100); avail: 32'
    }

    def 'should not fail when a task is within the total'() {
        when:
        new ResourceTracker(8, 0).validate(handler(8, '1TB'))
        then:
        noExceptionThrown()
    }

    def 'should sum the resources of a job array'() {
        given:
        def children = (1..3).collect { handler(2, '2GB') }
        def array = Mock(TaskHandler) {
            getTask() >> Mock(TaskArrayRun) { getChildren() >> children }
        }

        expect:
        ResourceTracker.cpus(array) == 6
        ResourceTracker.memory(array) == MemoryUnit.of('6GB').toBytes()
    }

    def 'should fail with negative values'() {
        when:
        new ResourceTracker(-1, 0)
        then:
        thrown(AssertionError)

        when:
        new ResourceTracker(0, -1)
        then:
        thrown(AssertionError)
    }
}
