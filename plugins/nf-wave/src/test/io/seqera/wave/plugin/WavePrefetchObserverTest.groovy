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

import groovyx.gpars.dataflow.DataflowQueue
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
}
