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

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
/**
 * Render the DAG using the Mermaid format to the specified file.
 * See https://mermaid-js.github.io/mermaid/#/ for more info.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class MermaidRenderer implements DagRenderer {

    static private final String INPUTS = 'inputs'

    static private final String OUTPUTS = 'outputs'

    private Session session = Global.session as Session

    private int depth = session.config.navigate('dag.depth', -1) as int

    private String direction = session.config.navigate('dag.direction', 'TB')

    private boolean verbose = session.config.navigate('dag.verbose', false);

    {
        if( direction !in ['TB','LR'] ) {
            log.warn "Invalid configuration property `dag.direction = '$direction'` - use either: 'TB' (top-bottom) or 'LR' (left-right)"
            this.direction = 'TB'
        }
    }

    @Override
    void renderDocument(DAG dag, Path file) {
        file.text = renderNetwork(dag)
    }

    String renderNetwork(DAG dag) {
        // construct node lookup from DAG
        def nodeLookup = getNodeLookup(dag)

        // collapse operator nodes
        if( !verbose )
            collapseOperators(nodeLookup)

        // remove empty workflow inputs
        if( !verbose )
            removeEmptyInputs(nodeLookup)

        // construct node tree
        final nodeTree = getNodeTree(nodeLookup)

        // collapse node tree to desired depth
        if( depth >= 0 )
            collapseNodeTree(nodeTree, depth, nodeLookup)

        // render diagram
        def lines = [] as List<String>
        lines << "flowchart ${direction}".toString()

        // render nodes
        renderNodeTree(lines, null, nodeTree)

        // render edges
        def edges = [] as Set<Edge>

        for( def node : nodeLookup.values() ) {
            for( def source : node.inputs )
                edges << new Edge(source.vertex, node.vertex)
            for( def target : node.outputs )
                edges << new Edge(node.vertex, target.vertex)
        }

        for( def edge : edges ) {
            final source = edge.source
            final target = edge.target
            def label = ""

            if( verbose ) {
                final dagEdges = dag.edges.findAll( e -> e.from == source && e.to == target && e.label )
                if( dagEdges )
                    label = "|${dagEdges*.label.join(',')}|".toString()
            }

            lines << "    ${source.name} -->${label} ${target.name}".toString()
        }

        lines << ""

        return lines.join('\n')
    }

    /**
     * Construct a map of nodes with inputs and outputs.
     *
     * @param dag
     */
    private Map<DAG.Vertex,Node> getNodeLookup(DAG dag) {
        def nodeLookup = [:] as Map<DAG.Vertex,Node>

        // create a node for each DAG vertex
        for( def v : dag.vertices ) {
            nodeLookup[v] = new Node(v)
        }

        // save the inputs and outputs for each node
        for( def e : dag.edges ) {
            final source = nodeLookup[e.from]
            final target = nodeLookup[e.to]
            source.outputs << target
            target.inputs << source

            // save labels for workflow inputs and outputs
            if( !verbose && e.from.type == DAG.Type.ORIGIN )
                source.label = e.label
            if( !verbose && e.to.type == DAG.Type.NODE )
                target.label = e.label
        }

        return nodeLookup
    }

    /**
     * Collapse operator nodes by replacing each operator subgraph
     * with a summary node.
     *
     * @param nodeLookup
     */
    private void collapseOperators(Map<DAG.Vertex,Node> nodeLookup) {
        def queue = nodeLookup
                .values()
                .findAll( n -> n.vertex.type == DAG.Type.OPERATOR ) as List<Node>

        while( !queue.isEmpty() ) {
            final node = queue.pop()
            final subgraph = findSubgraph(node, (Node n) -> n.vertex.type == DAG.Type.OPERATOR)
            collapseSubgraph(nodeLookup, subgraph, node.vertex)
            queue.removeAll(subgraph.nodes)
        }
    }

    /**
     * Remove 'Channel.empty' workflow inputs.
     *
     * @param nodeLookup
     */
    private void removeEmptyInputs(Map<DAG.Vertex,Node> nodeLookup) {
        def queue = nodeLookup
                .values()
                .findAll( n -> n.vertex.type == DAG.Type.ORIGIN )
                .findAll( n -> n.vertex.label == 'Channel.empty' ) as List<Node>

        while( !queue.isEmpty() ) {
            final node = queue.pop()
            final subgraph = findSubgraph(node, (Node n) -> n.vertex.type != DAG.Type.PROCESS)

            final hasProcessNeighbors = (subgraph.inputs + subgraph.outputs).any( (Node n) -> n.vertex.type == DAG.Type.PROCESS )
            if( hasProcessNeighbors )
                removeSubgraph(nodeLookup, new Subgraph([node], node.inputs, node.outputs))
            else {
                removeSubgraph(nodeLookup, subgraph)
                queue.removeAll(subgraph.nodes)
            }
        }
    }

    /**
     * Find the subgraph of all nodes connected to a seed node
     * that satisfy a predicate.
     *
     * @param seedNode
     * @param predicate
     */
    private Subgraph findSubgraph(Node seedNode, Closure predicate) {
        def queue = [seedNode] as List<Node>
        def subgraph = new Subgraph()

        while( !queue.isEmpty() ) {
            final v = queue.pop()

            for( def w : v.inputs + v.outputs ) {
                if( w in subgraph.nodes )
                    continue

                if( predicate(w) )
                    queue << w
                else if( w in v.inputs )
                    subgraph.inputs << w
                else  // w in v.outputs
                    subgraph.outputs << w
            }

            subgraph.nodes << v
        }

        return subgraph
    }

    /**
     * Collapse a subgraph by replacing it with a summary node.
     *
     * @param nodeLookup
     * @param subgraph
     * @param vertex
     */
    private void collapseSubgraph(Map<DAG.Vertex,Node> nodeLookup, Subgraph subgraph, DAG.Vertex vertex) {
        // remove subgraph
        removeSubgraph(nodeLookup, subgraph)

        // add summary node
        final summaryNode = new Node(vertex, null, subgraph.inputs, subgraph.outputs)
        nodeLookup[vertex] = summaryNode

        for( def w : subgraph.inputs )
            w.outputs << summaryNode

        for( def w : subgraph.outputs )
            w.inputs << summaryNode
    }

    /**
     * Remove a subgraph of nodes.
     *
     * @param nodeLookup
     * @param subgraph
     */
    private void removeSubgraph(Map<DAG.Vertex,Node> nodeLookup, Subgraph subgraph) {
        for( def w : subgraph.nodes )
            nodeLookup.remove(w.vertex)

        for( def w : subgraph.inputs )
            w.outputs -= subgraph.nodes

        for( def w : subgraph.outputs )
            w.inputs -= subgraph.nodes
    }

    /**
     * Construct a node tree with a subgraph for each subworkflow.
     *
     * @param nodeLookup
     */
    private Map<String,Object> getNodeTree(Map<DAG.Vertex,Node> nodeLookup) {
        // infer subgraphs of operator nodes
        final inferredKeys = inferSubgraphKeys(nodeLookup)

        // construct node tree
        def nodeTree = [:] as Map<String,Object>

        for( def node : nodeLookup.values() ) {
            final vertex = node.vertex

            // determine the vertex subgraph
            def keys = [] as List<String>

            if( vertex.type == DAG.Type.PROCESS ) {
                // extract keys from fully qualified name
                final result = getSubgraphKeys(vertex.label)
                keys = (List<String>)result[0]
                node.label = (String)result[1]
            }
            else if( vertex.type == DAG.Type.OPERATOR ) {
                // use inferred subgraph keys
                keys = inferredKeys[node]
            }
            else if( vertex.type == DAG.Type.ORIGIN ) {
                keys = [INPUTS]
            }
            else if( vertex.type == DAG.Type.NODE ) {
                keys = [OUTPUTS]
            }

            // navigate to given subgraph
            def subgraph = nodeTree
            for( def key : keys ) {
                if( key !in subgraph )
                    subgraph[key] = [:]
                subgraph = subgraph[key]
            }

            // add vertex to tree
            subgraph[vertex.name] = node
        }

        return nodeTree
    }

    /**
     * Infer the subgraph of each operator node from neighboring
     * processes.
     *
     * @param nodeLookup
     */
    private Map<Node,List> inferSubgraphKeys(Map<DAG.Vertex,Node> nodeLookup) {
        def inferredKeys = [:] as Map<Node,List>
        def queue = nodeLookup
                .values()
                .findAll( n -> n.vertex.type == DAG.Type.OPERATOR ) as List<Node>

        while( !queue.isEmpty() ) {
            // find subgraph of operator nodes
            final node = queue.pop()
            final subgraph = findSubgraph(node, (Node n) -> n.vertex.type == DAG.Type.OPERATOR)

            // select a neighboring process
            final inputs = subgraph.inputs.findAll( n -> n.vertex.type == DAG.Type.PROCESS )
            final outputs = subgraph.outputs.findAll( n -> n.vertex.type == DAG.Type.PROCESS )

            Node process = null
            if( inputs.size() == 1 )
                process = inputs[0]
            else if( outputs.size() == 1 )
                process = outputs[0]
            else if( inputs.size() > 0 )
                process = inputs[0]
            else if( outputs.size() > 0 )
                process = outputs[0]

            // extract keys from fully qualified process name
            final keys = process
                ? getSubgraphKeys(process.vertex.label)[0] as List
                : []

            // save inferred keys
            for( def w : subgraph.nodes )
                inferredKeys[w] = keys

            // update queue
            queue.removeAll(subgraph.nodes)
        }

        return inferredKeys
    }

    /**
     * Get the subgraph keys from a fully qualified process name.
     *
     * @param name
     */
    private List getSubgraphKeys(String name) {
        final tokens = name.tokenize(':')
        return [
            tokens[0..<tokens.size()-1],
            tokens.last()
        ]
    }

    /**
     * Collapse a node tree to the desired depth.
     *
     * @param nodeTree
     * @param depth
     * @param nodeLookup
     */
    private void collapseNodeTree(Map<String,Object> nodeTree, int depth, Map<DAG.Vertex,Node> nodeLookup) {
        nodeTree.each { key, value ->
            if( value !instanceof Map || key in [INPUTS, OUTPUTS] )
                return

            if( depth > 0 ) {
                collapseNodeTree((Map)value, depth - 1, nodeLookup)
                return
            }

            // collect subgraph
            final nodes = collectSubtree((Map)value)
            final inputs = nodes
                    .collect( n -> n.inputs )
                    .flatten()
                    .findAll( n -> n !in nodes ) as Set<Node>
            final outputs = nodes
                    .collect( n -> n.outputs )
                    .flatten()
                    .findAll( n -> n !in nodes ) as Set<Node>

            // remove subgraph
            removeSubgraph(nodeLookup, new Subgraph(nodes, inputs, outputs))

            // add summary node
            final vertex = nodes.first().vertex
            final summaryNode = new Node(vertex, key, inputs, outputs)
            nodeTree[key] = summaryNode
            nodeLookup[vertex] = summaryNode

            for( def w : inputs )
                w.outputs << summaryNode

            for( def w : outputs )
                w.inputs << summaryNode
        }
    }

    /**
     * Collect the nodes in a subtree.
     *
     * @param nodeTree
     */
    private List<Node> collectSubtree(Map<String,Object> nodeTree) {
        nodeTree.collect { key, value ->
                    value instanceof Map
                        ? collectSubtree(value)
                        : value
                }
                .flatten() as List<Node>
    }

    /**
     * Render a tree of nodes and subgraphs.
     *
     * @param lines
     * @param name
     * @param nodeTree
     */
    private void renderNodeTree(List<String> lines, String name, Map<String,Object> nodeTree) {
        if( name ) {
            final label = name in [INPUTS, OUTPUTS] ? '" "' : name
            lines << "    subgraph ${label}".toString()
        }

        nodeTree.each { key, value ->
            if( value instanceof Map )
                renderNodeTree(lines, key, value)
            else if( value instanceof Node ) {
                // skip node if it is disconnected
                if( !value.inputs && !value.outputs )
                    return
                lines << "    ${renderNode(value)}".toString()
            }
        }

        if( name )
            lines << "    end"
    }

    /**
     * Render a node.
     *
     * @param node
     */
    private String renderNode(Node node) {
        final id = node.vertex.name

        switch( node.vertex.type ) {
            case DAG.Type.PROCESS:
                return "${id}([${node.label}])"

            case DAG.Type.OPERATOR:
                return verbose
                    ? "${id}([${node.vertex.label}])"
                    : "${id}(( ))"

            case DAG.Type.ORIGIN:
            case DAG.Type.NODE:
                final label = node.vertex.label ?: node.label ?: ' '
                return "${id}[\"${label}\"]"

            default:
                return null
        }
    }

    @TupleConstructor
    private static class Node {
        DAG.Vertex vertex
        String label
        Set<Node> inputs = []
        Set<Node> outputs = []
    }

    @EqualsAndHashCode
    @TupleConstructor
    private static class Edge {
        DAG.Vertex source
        DAG.Vertex target
    }

    @TupleConstructor
    private static class Subgraph {
        List<Node> nodes = []
        Set<Node> inputs = []
        Set<Node> outputs = []
    }
}
