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

import nextflow.util.MemoryUnit
import spock.lang.Specification

class ResourceAccountTest extends Specification {

    private TaskHandler handler(int cpus, String memory) {
        Mock(TaskHandler) {
            getTask() >> new TaskRun(config: new TaskConfig(cpus: cpus, memory: MemoryUnit.of(memory)))
        }
    }

    def 'should be unlimited when created with zero values'() {
        when:
        def account = new ResourceAccount()
        then:
        account.isUnlimited()
        account.canReserve(100, MemoryUnit.of('100GB').toBytes())
        account.canReserve(handler(100, '100GB'))
    }

    def 'should track reserved resources'() {
        given:
        def account = new ResourceAccount(8, MemoryUnit.of('16GB').toBytes())
        and:
        def h = handler(4, '8GB')

        when:
        def reservation = account.reserve(h)
        then:
        !account.isUnlimited()
        account.availableCpus() == 4
        account.availableMemory() == MemoryUnit.of('8GB').toBytes()

        when: 'a task requiring the remaining resources is requested'
        then:
        account.canReserve(h)

        when: 'the remaining resources are reserved'
        account.reserve(h)
        then:
        account.availableCpus() == 0
        account.availableMemory() == 0
        !account.canReserve(h)
        !account.canReserve(handler(1, '1GB'))

        when: 'a reservation is released'
        account.release(reservation)
        then:
        account.availableCpus() == 4
        account.availableMemory() == MemoryUnit.of('8GB').toBytes()
        account.canReserve(h)
    }

    def 'should not allow a request exceeding the total'() {
        given:
        def account = new ResourceAccount(8, MemoryUnit.of('16GB').toBytes())

        expect:
        !account.canReserve(9, 0)
        !account.canReserve(0, MemoryUnit.of('17GB').toBytes())
    }

    def 'should support a limit on cpus only'() {
        given:
        def account = new ResourceAccount(4, 0)
        expect:
        account.canReserve(2, MemoryUnit.of('1TB').toBytes())
        !account.canReserve(5, 0)
    }

    def 'should support a limit on memory only'() {
        given:
        def account = new ResourceAccount(0, MemoryUnit.of('16GB').toBytes())
        expect:
        account.canReserve(100, MemoryUnit.of('8GB').toBytes())
        !account.canReserve(0, MemoryUnit.of('17GB').toBytes())
    }

    def 'should sum resources of a job array'() {
        given:
        def children = (1..3).collect { handler(2, '2GB') }
        def array = Mock(TaskHandler) {
            getTask() >> Mock(TaskArrayRun) { getChildren() >> children }
        }

        expect:
        ResourceAccount.cpusOf(array) == 6
        ResourceAccount.memOf(array) == MemoryUnit.of('6GB').toBytes()
    }

    def 'should fail with negative values'() {
        when:
        new ResourceAccount(-1, 0)
        then:
        thrown(AssertionError)

        when:
        new ResourceAccount(0, -1)
        then:
        thrown(AssertionError)
    }
}
