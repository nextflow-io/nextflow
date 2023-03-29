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
import spock.lang.Specification

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class MermaidRendererTest extends Specification {

    def setupSpec() {
        new Session()
    }

    def 'should render an abstract graph using the `mmd` format' () {
        given:
        def file = Files.createTempFile('test', null)
        def ch1 = new DataflowQueue()
        def ch2 = new DataflowQueue()
        def ch3 = new DataflowQueue()

        def dag = new DAG()
        dag.addOperatorNode('Op1', ch1, ch2)
        dag.addOperatorNode('Op2', ch2, ch3)

        dag.normalize()

        when:
        new MermaidRenderer().renderAbstractGraph(dag, file)
        then:
        file.text ==
            '''
            flowchart TD
                p0(( ))
                p1([Op1])
                p2([Op2])
                p3(( ))
                p0 --> p1
                p1 --> p2
                p2 --> p3
            '''
            .stripIndent().leftTrim()

        cleanup:
        file.delete()
    }

    def 'should render a concrete graph using the `mmd` format' () {
        given:
        def file = Files.createTempFile('test', null)

        def dag = Mock(ConcreteDAG) {
            nodes >> [
                '012345': new ConcreteDAG.Task(
                    index: 1,
                    label: 'foo',
                    inputs: [
                        new ConcreteDAG.Input(
                            name: 'data.txt',
                            path: Paths.get('/inputs/data.txt'),
                            predecessor: null
                        )
                    ],
                    outputs: [
                        new ConcreteDAG.Output(
                            name: 'data.foo',
                            path: Paths.get('/work/012345/data.foo'),
                        )
                    ]
                ),
                'abcdef': new ConcreteDAG.Task(
                    index: 2,
                    label: 'bar',
                    inputs: [
                        new ConcreteDAG.Input(
                            name: 'data.foo',
                            path: Paths.get('/work/012345/data.foo'),
                            predecessor: '012345'
                        )
                    ],
                    outputs: [
                        new ConcreteDAG.Output(
                            name: 'data.bar',
                            path: Paths.get('/work/abcdef/data.bar'),
                        )
                    ]
                )
            ]
        }

        when:
        new MermaidRenderer().renderConcreteGraph(dag, file)
        then:
        file.text ==
            '''
            flowchart TD
                t1["foo"]
                i1(( )) -->|data.txt| t1
                t2["bar"]
                t1 -->|data.foo| t2
                t2 -->|data.bar| o1(( ))
            '''
            .stripIndent().leftTrim()

        cleanup:
        file.delete()
    }
}
