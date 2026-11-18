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
package nextflow.dataflow.engine

import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import spock.lang.Specification

/**
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ProcessBacklogTest extends Specification {

    private static ProcessBacklog backlog(Backpressure bp, DataflowReadChannel... inputs) {
        return new ProcessBacklog(bp, 'foo', inputs.toList(), null)
    }

    def 'should close the operator gate while preparing and pending tasks reach the high-water mark' () {
        given:
        def backlog = backlog(new Backpressure(null, 8))

        when:
        3.times { backlog.taskStarted() }
        5.times { backlog.taskPending() }
        then:
        backlog.getPreparing() == 3
        backlog.getPending() == 5
        backlog.getAsBoolean()

        when:
        backlog.taskPrepared()
        then:
        !backlog.getAsBoolean()

        when:
        backlog.taskPending()
        backlog.close()
        then:
        !backlog.getAsBoolean()
    }

    def 'should open the operator gate when backpressure is released' () {
        given:
        def bp = new Backpressure(null, 2)
        def backlog = backlog(bp)

        when:
        2.times { backlog.taskPending() }
        then:
        backlog.getAsBoolean()

        when:
        bp.release()
        then:
        !backlog.getAsBoolean()
    }

    def 'should count the complete input sets' () {
        given:
        def a = new DataflowBroadcast()
        def b = new DataflowBroadcast()
        def q = new DataflowQueue()
        def backlog = backlog(new Backpressure(null, 8), a.createReadChannel(), b.createReadChannel(), q)
        assert backlog.ready(100) == 0

        when:
        for( int i = 0; i < 5; i++ ) a.bind(i)
        for( int i = 0; i < 2; i++ ) b.bind(i)
        for( int i = 0; i < 3; i++ ) q.bind(i)
        then:
        backlog.ready(100) == 2
        backlog.ready(1) == 1
    }

    def 'should not count items waiting for an unbound value' () {
        given:
        def ch = new DataflowBroadcast()
        def value = new DataflowVariable()
        def backlog = backlog(new Backpressure(null, 4), ch.createReadChannel(), value)

        when:
        for( int i = 0; i < 10; i++ ) ch.bind(i)
        then:
        backlog.ready(100) == 0
        !backlog.isOver(4)

        when:
        value.bind('x')
        then:
        backlog.ready(100) == 10
        backlog.isOver(4)
    }

    def 'should not count a process with only value inputs' () {
        given:
        def value = new DataflowVariable()
        value.bind(1)
        def backlog = backlog(new Backpressure(null, 4), value)

        expect:
        backlog.ready(100) == 0
    }

    def 'should compare the whole backlog with a mark' () {
        given:
        def ch = new DataflowBroadcast()
        def backlog = backlog(new Backpressure(null, 8), ch.createReadChannel())

        when:
        2.times { backlog.taskStarted() }
        3.times { backlog.taskPending() }
        ch.bind(0)
        ch.bind(1)
        then:
        !backlog.isOver(8)
        backlog.isOver(7)

        when:
        backlog.close()
        then:
        !backlog.isOver(1)
    }
}
