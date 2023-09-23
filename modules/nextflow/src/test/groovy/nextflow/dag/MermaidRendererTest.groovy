/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.dag

import java.nio.file.Files
import java.nio.file.Paths

import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Session
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class MermaidRendererTest extends Specification {

    def setupSpec() {
        new Session()
    }

    def 'should render a workflow graph using the `mmd` format' () {
        given:
        def file = Files.createTempFile('test', null)
        def ch1 = new DataflowQueue()
        def ch2 = new DataflowQueue()
        def ch3 = new DataflowQueue()

        def session = new Session([dag: [verbose: true]])
        def dag = new DAG()
        dag.addOperatorNode('Op1', ch1, ch2)
        dag.addOperatorNode('Op2', ch2, ch3)

        dag.normalize()

        when:
        new MermaidRenderer().renderWorkflowGraph(dag, file)
        then:
        file.text ==
            '''
            flowchart TD
                subgraph " "
                v0[" "]
                end
                v1([Op1])
                v2([Op2])
                subgraph " "
                v3[" "]
                end
                v0 --> v1
                v1 --> v2
                v2 --> v3
            '''
            .stripIndent().trim()

        cleanup:
        file.delete()
    }

    def 'should render a task graph using the `mmd` format' () {
        given:
        def file = Files.createTempFile('test', null)

        def task1 = Mock(TaskRun)
        def task2 = Mock(TaskRun)
        def output1 = Paths.get('/work/012345/data.foo')
        def v1 = new TaskDAG.Vertex(
            index: 1,
            label: 'foo',
            inputs: [ 'data.txt': Paths.get('/inputs/data.txt') ],
            outputs: [ output1 ]
        )
        def v2 = new TaskDAG.Vertex(
            index: 2,
            label: 'bar',
            inputs: [ 'data.foo': output1 ],
            outputs: [ Paths.get('/work/abcdef/data.bar') ]
        )
        def dag = Mock(TaskDAG) {
            vertices >> [
                (task1): v1,
                (task2): v2
            ]
            getProducerVertex(output1) >> v1
        }

        when:
        new MermaidRenderer().renderTaskGraph(dag, file)
        then:
        file.text ==
            '''
            flowchart TD
                t1["foo"]
                i1(( )) -->|/inputs/data.txt| t1
                t2["bar"]
                t1 -->|data.foo| t2
                t2 -->|data.bar| o1(( ))
            '''
            .stripIndent().trim()

        cleanup:
        file.delete()
    }
}
