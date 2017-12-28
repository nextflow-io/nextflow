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
import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Global
import nextflow.Session
import nextflow.processor.TaskProcessor
import nextflow.script.InputsList
import nextflow.script.OutputsList
/**
 * Helper class to mark DAG node with the proper labels
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class NodeMarker {

    static private List<DataflowProcessor> operators = []

    static void appendOperator( DataflowProcessor p ) {
        operators.add(p)
    }

    static private Session getSession() { Global.session as Session }

    /**
     *  Creates a new vertex in the DAG representing a computing `process`
     *
     * @param label The label associated to the process
     * @param inputs The list of inputs entering in the process
     * @param outputs the list of outputs leaving the process
     */
    static void addProcessNode( TaskProcessor process, InputsList inputs, OutputsList outputs ) {
        if( session && session.dag && !session.aborted )
            session.dag.addProcessNode( process.name, inputs, outputs, process )
    }

    /**
     * Creates a new DAG vertex representing a dataflow operator
     *
     * @param label The operator label
     * @param inputs The operator input(s). It can be either a single channel or a list of channels.
     * @param outputs The operator output(s). It can be either a single channel, a list of channels or {@code null} if the operator has no output.
     */
    static void addOperatorNode( String name, inputs, outputs )  {
        if( session && session.dag && !session.aborted )
            session.dag.addOperatorNode(name, inputs, outputs, operators)
    }

    /**
     * Creates a vertex in the DAG representing a dataflow channel source.
     *
     * @param label The node description
     * @param source Either a dataflow channel or a list of channel.
     */
    static void addSourceNode( String name, source )  {
        if( session && session.dag && !session.aborted )
            session.dag.addSourceNode(name, source)
    }


}
