/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
            digraph TheGraph {
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
