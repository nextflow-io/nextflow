/*
 * Copyright 2013-2024, Seqera Labs
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

import nextflow.NextflowMeta
import nextflow.processor.TaskProcessor
import nextflow.script.BaseScript
import nextflow.script.ProcessConfig
import nextflow.script.params.InParam
import nextflow.script.params.InputsList
import nextflow.script.params.OutParam
import nextflow.script.params.OutputsList

import java.nio.file.Files

import groovyx.gpars.dataflow.DataflowQueue
import nextflow.trace.config.DagConfig
import spock.lang.Specification
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class MermaidRendererTest extends Specification {

    //This test requires DSL2 to be tested alone
    void setup(){
        NextflowMeta.instance.enableDsl2()
    }

    def 'should render a graph using the `mmd` format' () {
        given:
        def file = Files.createTempFile('test', null)
        def ch1 = new DataflowQueue()
        def ch2 = new DataflowQueue()
        def ch3 = new DataflowQueue()
        and:
        def dag = new DAG()
        dag.addOperatorNode('Op1', ch1, ch2)
        dag.addOperatorNode('Op2', ch2, ch3)
        dag.normalize()
        and:
        def config = new DagConfig(verbose: true)

        when:
        new MermaidRenderer(config).renderDocument(dag, file)
        then:
        file.text ==
            '''
            flowchart TB
                subgraph " "
                v0[" "]
                end
                v1(["Op1"])
                v2(["Op2"])
                subgraph " "
                v3[" "]
                end
                v0 --> v1
                v1 --> v2
                v2 --> v3
            '''
            .stripIndent().leftTrim()

        cleanup:
        file.delete()
    }

    def 'should render a graph with description using the `mmd` format' () {
        given:
        def file = Files.createTempFile('test', null)
        def ch1 = new DataflowQueue()
        def ch2 = new DataflowQueue()
        def ch3 = new DataflowQueue()

        def pInList = new InputsList()
        def ip1 = Mock(InParam)
        pInList.add( ip1 )

        def pOutList = new OutputsList()
        def op1 = Mock(OutParam)
        pOutList.add( op1 )

        def process = Mock(TaskProcessor){

            getConfig() >> new ProcessConfig(Mock(BaseScript), "test")
        }
        process.config.meta([description: 'a short description of the process'])

        and:
        def dag = new DAG()
        dag.addOperatorNode('Op1', ch1, ch2)
        dag.addOperatorNode('Op2', ch2, ch3)
        dag.addProcessNode('a process', pInList, pOutList, process)
        dag.normalize()
        and:
        def config = new DagConfig(verbose: true)

        when:
        new MermaidRenderer(config).renderDocument(dag, file)
        then:
        file.text ==
            '''
            flowchart TB
                subgraph " "
                v0[" "]
                v4[" "]
                end
                v1(["Op1"])
                v2(["Op2"])
                subgraph " "
                v3[" "]
                end
                v5(["a process"])
                note_v5[a short description of the process]
                v5 --- note_v5
                style note_v5 fill:#ffffff,stroke:#000,stroke-width:3px,color:#000,stroke-dasharray: 5 5
                v0 --> v1
                v1 --> v2
                v2 --> v3
                v4 --> v5
            '''
                .stripIndent().leftTrim()

        cleanup:
        file.delete()
    }
}
