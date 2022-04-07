/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import groovyx.gpars.dataflow.DataflowQueue
import spock.lang.Specification

import nextflow.Session

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class MermaidRendererTest extends Specification {

    def setupSpec() {
        new Session()
    }

    def 'should render a graph using the `mmd` format' () {
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
        new MermaidRenderer().renderDocument(dag, file)
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
}
