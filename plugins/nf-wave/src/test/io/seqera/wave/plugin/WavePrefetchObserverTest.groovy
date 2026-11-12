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

package io.seqera.wave.plugin

import java.time.OffsetDateTime

import groovyx.gpars.dataflow.DataflowQueue
import nextflow.container.resolver.ContainerMeta
import nextflow.dag.DAG
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.params.InParam
import nextflow.script.params.InputsList
import nextflow.script.params.OutParam
import nextflow.script.params.OutputsList
import nextflow.trace.event.TaskEvent
import spock.lang.Specification
import spock.util.concurrent.PollingConditions

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WavePrefetchObserverTest extends Specification {

    static class RecordingObserver extends WavePrefetchObserver {
        List<String> submitted = []
        RecordingObserver(DAG dag, int depth) { super(dag, depth) }
        @Override protected void submit(TaskProcessor process) { submitted << process.name }
    }

    TaskProcessor process(String name) {
        Mock(TaskProcessor) { getName() >> name }
    }

    InputsList inputs(DataflowQueue... channels) {
        final result = new InputsList()
        for( def c : channels ) { final ch = c; result.add(Mock(InParam) { getRawChannel() >> ch }) }
        return result
    }

    OutputsList outputs(DataflowQueue... channels) {
        final result = new OutputsList()
        for( def c : channels ) { final ch = c; result.add(Mock(OutParam) { getOutChannel() >> ch }) }
        return result
    }

    TaskEvent pending(TaskProcessor process) {
        final task = new TaskRun(processor: process)
        new TaskEvent(Mock(TaskHandler) { getTask() >> task }, null)
    }

    /*
     *  A --> map --> B --> C
     *  A --> D
     */
    def dag(Map<String,TaskProcessor> p) {
        final ab1 = new DataflowQueue(); final ab2 = new DataflowQueue(); final bc = new DataflowQueue(); final ad = new DataflowQueue()
        final dag = new DAG()
        dag.addProcessNode('A', inputs(new DataflowQueue()), outputs(ab1, ad), p.A)
        dag.addOperatorNode('map', ab1, ab2)
        dag.addProcessNode('B', inputs(ab2), outputs(bc), p.B)
        dag.addProcessNode('C', inputs(bc), outputs(new DataflowQueue()), p.C)
        dag.addProcessNode('D', inputs(ad), outputs(new DataflowQueue()), p.D)
        return dag
    }

    def 'should find the downstream processes up to the given depth' () {
        given:
        def p = [A: process('A'), B: process('B'), C: process('C'), D: process('D')]
        def dag = dag(p)

        expect:
        WavePrefetchObserver.downstream(dag, p[FROM], DEPTH)*.name as Set == EXPECTED as Set

        where:
        FROM | DEPTH | EXPECTED
        'A'  | 1     | ['B', 'D']
        'A'  | 2     | ['B', 'D', 'C']
        'B'  | 2     | ['C']
        'C'  | 2     | []
    }

    def 'should request each downstream container once' () {
        given:
        def p = [A: process('A'), B: process('B'), C: process('C'), D: process('D')]
        def observer = new RecordingObserver(dag(p), 1)

        when: 'the first task of A is pending'
        observer.onTaskPending(pending(p.A))
        then:
        observer.submitted as Set == ['B', 'D'] as Set

        when: 'more tasks of A are pending'
        observer.onTaskPending(pending(p.A))
        then:
        observer.submitted.size() == 2

        when: 'the first task of B is pending: C is requested, B is not requested again'
        observer.onTaskPending(pending(p.B))
        then:
        observer.submitted.size() == 3
        observer.submitted.last() == 'C'
    }

    def 'should not request a process whose own task came first' () {
        given:
        def p = [A: process('A'), B: process('B'), C: process('C'), D: process('D')]
        def observer = new RecordingObserver(dag(p), 1)

        when:
        observer.onTaskPending(pending(p.D))
        observer.onTaskPending(pending(p.A))
        then:
        observer.submitted == ['B']
    }

    def 'should request a container already ready' () {
        given:
        def observer = new WavePrefetchObserver(new DAG(), 1)
        def task = Mock(TaskRun)
        def B = Mock(TaskProcessor) { getName() >> 'B' }

        when:
        observer.prefetch(B)
        then:
        1 * B.createTaskPreview() >> task
        1 * task.getContainer() >> 'wave.io/wt/xyz/foo:latest'
        1 * task.isContainerReady() >> true
        and:
        observer.countReady.get() == 1
        observer.countBuilding.get() == 0
    }

    def 'should not await a container build' () {
        given:
        def observer = new WavePrefetchObserver(new DAG(), 1)
        def task = Mock(TaskRun)
        def B = Mock(TaskProcessor) { getName() >> 'B' }

        when:
        observer.prefetch(B)
        then: 'the container status is checked once'
        1 * B.createTaskPreview() >> task
        1 * task.getContainer() >> 'wave.io/wt/xyz/foo:latest'
        1 * task.isContainerReady() >> false
        and:
        observer.countBuilding.get() == 1
        observer.countReady.get() == 0
    }

    def 'should skip a queued process whose own task already started' () {
        given:
        def p = [A: process('A'), B: process('B'), C: process('C'), D: process('D')]
        def observer = new WavePrefetchObserver(dag(p), 1)

        when:
        observer.onTaskPending(pending(p.B))
        observer.prefetch(p.B)
        then:
        0 * p.B.createTaskPreview()
    }

    def 'should send the requests from the prefetch thread' () {
        given:
        def observer = new WavePrefetchObserver(new DAG(), 1)
        def threads = Collections.synchronizedList([])
        def total = 3
        def processes = (1..total).collect { int i ->
            final task = Stub(TaskRun) { isContainerReady() >> false }
            Stub(TaskProcessor) {
                getName() >> "P$i"
                createTaskPreview() >> { threads << WavePrefetchObserver.isPrefetchThread(); task }
            }
        }
        def conditions = new PollingConditions(timeout: 5)

        when: 'all the containers are building'
        processes.each { observer.submit(it) }
        then: 'all of them are requested without awaiting the builds'
        conditions.eventually { assert observer.countBuilding.get() == total }
        and:
        threads == [true]*total
        !WavePrefetchObserver.isPrefetchThread()

        cleanup:
        observer.onFlowComplete()
    }

    def 'should describe the container request' () {
        given:
        def now = OffsetDateTime.now()
        def meta = new ContainerMeta(requestId: '12345', buildId: 'bd-123_1', cached: true, scanId: 'sc-123_1', sourceImage: 'cr.io/wave/build:foo', requestTime: now)
        def task = Stub(TaskRun) { containerMeta() >> meta }
        def noMeta = Stub(TaskRun) { containerMeta() >> null }
        def startedAt = now.toInstant().toEpochMilli()

        expect: 'a response received after the prefetch started'
        WavePrefetchObserver.describe(task, startedAt) == 'response: new; requestId: 12345; buildId: bd-123_1; cached: true; scanId: sc-123_1; source: cr.io/wave/build:foo'
        and: 'a response received before is reused from the Wave client cache'
        WavePrefetchObserver.describe(task, startedAt + 1).startsWith('response: reused; requestId: 12345;')
        and:
        WavePrefetchObserver.describe(noMeta, startedAt) == 'meta: none'
    }
}
