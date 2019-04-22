/*
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
import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Global
import nextflow.Session
import nextflow.processor.TaskProcessor
import nextflow.script.params.InputsList
import nextflow.script.params.OutputsList
/**
 * Helper class to mark DAG node with the proper labels
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class NodeMarker {

    static private List<DataflowProcessor> operators = new ArrayList<>(10)

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
        if( session && session.dag && !session.aborted ) {
            session.dag.addOperatorNode(name, inputs, outputs, new ArrayList<DataflowProcessor>(operators))
            operators.clear()
        }
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
