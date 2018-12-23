/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DotRendererTest extends Specification {

    def 'should remove not alphanumeric chars' (){

        expect:
        DotRenderer.normalise('hello') == 'hello'
        DotRenderer.normalise('qwe ,r/t`y*$*$')  == 'qwerty'
    }

    def 'should render a graph y using the `dot` format' () {

        given:
        def file = Files.createTempFile('test',null)
        def ch1 = new DataflowQueue()
        def ch2 = new DataflowQueue()
        def ch3 = new DataflowQueue()

        def dag = new DAG()
        dag.addOperatorNode('Op1', ch1, ch2)
        dag.addOperatorNode('Op2', ch2, ch3)

        dag.normalize()

        when:
        new DotRenderer('TheGraph').renderDocument(dag, file)
        then:
        file.text ==
            '''
            digraph "TheGraph" {
            p0 [shape=point,label="",fixedsize=true,width=0.1];
            p1 [shape=circle,label="",fixedsize=true,width=0.1,xlabel="Op1"];
            p0 -> p1;

            p1 [shape=circle,label="",fixedsize=true,width=0.1,xlabel="Op1"];
            p2 [shape=circle,label="",fixedsize=true,width=0.1,xlabel="Op2"];
            p1 -> p2;

            p2 [shape=circle,label="",fixedsize=true,width=0.1,xlabel="Op2"];
            p3 [shape=point];
            p2 -> p3;

            }
            '''
            .stripIndent().leftTrim()

        cleanup:
        file.delete()


    }
}
