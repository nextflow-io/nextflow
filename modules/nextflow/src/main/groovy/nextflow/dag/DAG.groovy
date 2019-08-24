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

import groovy.transform.MapConstructor
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.NF
import nextflow.extension.CH
import nextflow.extension.DataflowHelper
import nextflow.processor.TaskProcessor
import nextflow.script.params.DefaultInParam
import nextflow.script.params.DefaultOutParam
import nextflow.script.params.EachInParam
import nextflow.script.params.InParam
import nextflow.script.params.InputsList
import nextflow.script.params.OutParam
import nextflow.script.params.OutputsList
import nextflow.script.params.TupleInParam
import nextflow.script.params.TupleOutParam
/**
 * Model a direct acyclic graph of the pipeline execution.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class DAG {

    @PackageScope
    static enum Type {
        PROCESS,
        OPERATOR,
        ORIGIN,
        NODE
    }

    /**
     * The list of edges in the graph
     */
    private List<Edge> edges = new ArrayList<>(50)

    /**
     * The ordered list of vertices
     */
    private List<Vertex> vertices = new ArrayList<>(50)

    @PackageScope
    List<Vertex> getVertices() { vertices }

    @PackageScope
    List<Edge> getEdges() { edges }

    boolean isEmpty() { edges.size()==0 && vertices.size()==0 }

    /**
     *  Creates a new vertex in the DAG representing a computing `process`
     *
     * @param label The label associated to the process
     * @param inputs The list of inputs entering in the process
     * @param outputs the list of outputs leaving the process
     */
    void addProcessNode( String label, InputsList inputs, OutputsList outputs, TaskProcessor process=null ) {
        assert label
        assert inputs
        assert outputs
        addVertex( Type.PROCESS, label, normalizeInputs(inputs), normalizeOutputs(outputs), process )
    }

    /**
     * Creates a new DAG vertex representing a dataflow operator
     *
     * @param label The operator label
     * @param inputs The operator input(s). It can be either a single channel or a list of channels.
     * @param outputs The operator output(s). It can be either a single channel, a list of channels or {@code null} if the operator has no output.
     */
    void addOperatorNode( String label, inputs, outputs, List<DataflowProcessor> operators=null )  {
        assert label
        assert inputs
        addVertex(Type.OPERATOR, label, normalizeChannels(inputs), normalizeChannels(outputs), operators )
    }

    /**
     * Creates a vertex in the DAG representing a dataflow channel source.
     *
     * @param label The node description
     * @param source Either a dataflow channel or a list of channel.
     */
    void addSourceNode( String label, source )  {
        assert label
        assert source
        addVertex(Type.ORIGIN, label, null, normalizeChannels(source) )
    }

    /**
     * Creates a vertex and adds it to the DAG
     *
     * @param type A {link Type} value
     * @param label The vertex description
     * @param inbounds The inbounds channels to this vertex
     * @param outbounds The outbounds channels leaving the vertex
     */
    @PackageScope
    void addVertex( Type type, String label, List<ChannelHandler> inbounds, List<ChannelHandler> outbounds, Object extra=null) {

        final vertex = createVertex( type, label, extra )

        for( ChannelHandler channel : inbounds ) {
            inbound( vertex, channel )
        }

        for( ChannelHandler channel : outbounds ) {
            outbound( vertex, channel )
        }
    }

    /**
     * Creates a DAG vertex object
     *
     * @param type The vertex type
     * @param label The vertex label
     * @return A {@link Vertex} object
     */
    @PackageScope
    Vertex createVertex( Type type, String label, extra=null ) {
        def result = new Vertex(type, label)
        if( extra instanceof TaskProcessor ) {
            result.process = extra
            result.operators = [ extra.operator ]
        }
        else if( extra instanceof List ) {
            result.operators = (List)extra.clone()
        }
        else if( extra != null )
            throw new IllegalArgumentException("Not a valid DAG vertex parameter: [${extra.class.name}] $extra")

        vertices << result
        return result
    }

    private void inbound( Vertex vertex, ChannelHandler entering )  {

        // look for an existing edge for the given dataflow channel
        def edge = findEdge(entering.channel)

        // if does not exist just create it
        if( !edge ) {
            edges << new Edge(channel: entering.channel, to: vertex, label: entering.label)
        }
        // link the edge to given `edge`
        else if( edge.to == null ) {
            edge.to = vertex
        }
        // handle the special case for dataflow variable
        // this kind of channel can be used more than one time as an input
        else if( isForkable(entering.channel) ) {
            if( !edge.from ) {
                edge.from = new Vertex(Type.ORIGIN);
                int p = vertices.indexOf(edge.to)
                if(p!=-1) vertices.add(p,edge.from)
                else vertices.add(edge.from)
            }
            def fork = new Edge(channel: entering.channel, from: edge.from, to: vertex, label: entering.label)
            edges << fork
        }
        // the same channel - apart the above case - cannot be used multiple times as an input
        // thus throws an exception
        else {
            final name = getChannelName(entering)
            throw new MultipleInputChannelException(name, vertex, edge.to)
        }
    }

    private boolean isForkable(obj) {
        if( obj instanceof DataflowExpression )
            return true
        if( obj instanceof DataflowBroadcast )
            return true
        return obj instanceof DataflowQueue && CH.isBridge(obj)
    }

    private void outbound( Vertex vertex, ChannelHandler leaving) {

        // look for an existing edge for the given dataflow channel
        final edge = findEdge(leaving.channel)
        if( !edge ) {
            edges << new Edge(channel: leaving.channel, from: vertex, label: leaving.label)
        }
        else if( edge.from == null ) {
            edge.from = vertex
        }
        // the same channel cannot be used multiple times as an output
        // thus throws an exception
        else {
            final name = getChannelName(leaving)
            throw new MultipleOutputChannelException(name, vertex, edge.from)
        }

    }

    static private List<ChannelHandler> normalizeInputs( InputsList inputs ) {

        inputs
                .findAll { !( it instanceof DefaultInParam)  }
                .collect { InParam p -> new ChannelHandler(channel: p.rawChannel, label: inputName0(p)) }

    }

    static private String inputName0(InParam param) {
        if( param instanceof TupleInParam ) return null
        if( param instanceof EachInParam ) return null
        return param.name
    }

    static private List<ChannelHandler> normalizeOutputs( OutputsList outputs ) {

        def result = []
        for(OutParam p :outputs) {
            if( p instanceof DefaultOutParam ) break
            for(Object it : p.outChannels) {
                result << new ChannelHandler(channel: it, label: p instanceof TupleOutParam ? null : p.name)
            }
        }

        return result
    }

    static private List<ChannelHandler> normalizeChannels( entry ) {
        if( entry == null ) {
            Collections.emptyList()
        }
        else if( entry instanceof DataflowReadChannel || entry instanceof DataflowWriteChannel ) {
            [ new ChannelHandler(channel: entry) ]
        }
        else if( entry instanceof Collection || entry instanceof Object[] ) {
            entry.collect { new ChannelHandler(channel: it) }
        }
        else {
            throw new IllegalArgumentException("Not a valid channel type: [${entry.class.name}]")
        }
    }

    @PackageScope
    Edge findEdge( channel ) {
        edges.find { edge -> edge.channel.is(channel) }
    }

    @PackageScope
    int indexOf(Vertex v) {
        vertices.indexOf(v)
    }

    @PackageScope
    void normalizeMissingVertices() {
        for( Edge e : edges ) {
            assert e.from || e.to, 'Missing source and termination vertices for edge'

            if( !e.from ) {
                // creates the missing origin vertex
                def vertex = e.from = new Vertex(Type.ORIGIN)
                int p = vertices.indexOf( e.to )
                vertices.add( p, vertex )
            }
            else if( !e.to ) {
                // creates the missing termination vertex
                def vertex = e.to = new Vertex(Type.NODE)
                int p = vertices.indexOf( e.from )
                vertices.add( p+1, vertex )
            }
        }
    }

    @PackageScope
    void resolveEdgeNames() {
        for( Edge edge : edges ) {
            final name = lookupVariable(edge.channel)
            if( name )
                edge.label = name
        }
    }

    @PackageScope String lookupVariable(obj) {
        NF.lookupVariable(obj)
    }

    @PackageScope
    String resolveChannelName( Map map, channel ) {
        def entry = map.find { k,v -> v.is channel }
        return entry ? entry.key : null
    }

    @PackageScope
    String getChannelName( ChannelHandler handler ) {
        NF.lookupVariable(handler.channel) ?: handler.label
    }

    void normalize() {
        normalizeMissingVertices()
        resolveEdgeNames()
    }

    /**
     * @return
     *      A string listing the current active processes/operators in the
     *      dataflow network represented by this DAG
     */
    String dumpActiveNodes() {
        normalize()

        // first dump active processes
        def processes = vertices.findAll { it.process && it.isActive() }.collect { it.process }
        if( processes ) {
            def result = new StringBuilder()
            processes.eachWithIndex { it, index ->
                if( index>0 ) result << '\n'
                result << it.dumpTerminationStatus()
            }
            return result.toString()
        }

        // otherwise fallback on other nodes
        def nodes = vertices.findAll { it.active }
        if( !nodes )
            return null

        def result = new StringBuilder()
        nodes.each {
            result << '  [' << it.type.toString().toLowerCase() << "] " << (it.label ?: it.name) << '\n'
        }

        return result
    }

    /**
     * Model a vertex in the DAG.
     *
     * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
     */
    @ToString(includeNames = true, includes = 'label,type', includePackage=false)
    @PackageScope
    class Vertex {

        /**
         * The vertex label
         */
        String label

        /**
         * The vertex type
         */
        Type type

        /**
         * One or more {@link DataflowProcessor} associated to this graph node
         */
        List<DataflowProcessor> operators

        TaskProcessor process

        /**
         * Create an DGA vertex instance
         *
         * @param type A {@link Type} value
         * @param label A descriptive string to label this vertex
         */
        Vertex( Type type, String label = null ) {
            assert type
            this.label = label
            this.type = type
        }

        /**
         * @return The order of the index in the DAG
         */
        int getOrder() {
            indexOf(this)
        }

        /**
         * @return The unique name for this node
         */
        String getName() { "p${getOrder()}" }

        boolean isActive() {
            operators?.any { DataflowHelper.isProcessorActive(it) }
        }

    }

    /**
     * Models an edge in the DAG
     *
     * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
     */
    @PackageScope
    @ToString(includeNames = true, includes = 'label,from,to', includePackage=false)
    @MapConstructor
    class Edge {

        /**
         * The Dataflow channel that originated this graph edge
         */
        Object channel

        /**
         * The vertex *from* where the edge starts
         */
        Vertex from

        /**
         * The vertex *to* where the edge ends
         */
        Vertex to

        /**
         * A descriptive label
         */
        String label

    }

    /**
     * A simple wrapper object to handle a channel and the associated label
     */
    @ToString(includeNames = true, includes = 'label', includePackage=false)
    static class ChannelHandler {

        /**
         * The {@link groovyx.gpars.dataflow.DataflowChannel} that originated this graph edge
         */
        Object channel

        /**
         * The edge label
         */
        String label

    }


}
