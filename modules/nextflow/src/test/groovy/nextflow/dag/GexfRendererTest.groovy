/*
 * Copyright 2013-2020, Université de Nantes, CNRS, INSERM, l’institut du thorax, F-44000 Nantes, France.
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
import java.io.File

import groovyx.gpars.dataflow.DataflowQueue
import groovy.util.XmlSlurper
import spock.lang.Specification

import nextflow.Session

/**
 *
 * @author Pierre Lindenbaum / yokofakun
 */
class GexfRendererTest extends Specification {

    def setupSpec() {
        new Session()
    }

    def 'test graph is a xml document with -gexf- root ' () {

        given:
            def file = File.createTempFile("test",".xml")
            def ch1 = new DataflowQueue()
            def ch2 = new DataflowQueue()
            def ch3 = new DataflowQueue()

            def dag = new DAG()
            dag.addOperatorNode('Op1', ch1, ch2)
            dag.addOperatorNode('Op2', ch2, ch3)

            dag.normalize()
        when:
                new GexfRenderer('TheGraph').renderDocument(dag, file.toPath())
        then:
                def graph = new XmlSlurper().parse(file);
                assert graph.name() == 'gexf'
        cleanup:
                file.delete()

    }
}
