package nextflow.dag
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowChannel
import nextflow.Session
import nextflow.script.DefaultInParam
import nextflow.script.DefaultOutParam
import nextflow.script.InParam
import nextflow.script.InputsList
import nextflow.script.OutParam
import nextflow.script.OutputsList
import nextflow.script.SetInParam
import nextflow.script.SetOutParam
/**
 * Model the direct acyclic graph of a pipeline execution.
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

    /**
     * The {@link Session} to which this DAG is bound
     */
    private Session session

    @PackageScope
    List<Vertex> getVertices() { vertices }

    @PackageScope
    List<Edge> getEdges() { edges }

    /**
     *  Creates a new vertex in the DAG representing a computing `process`
     *
     * @param label The label associated to the process
     * @param inputs The list of inputs entering in the process
     * @param outputs the list of outputs leaving the process
     */
    void addProcessNode( String label, InputsList inputs, OutputsList outputs ) {
        assert label
        assert inputs
        assert outputs
        addVertex( Type.PROCESS, label, normalizeInputs(inputs), normalizeOutputs(outputs))

    }

    /**
     * Creates a new DAG vertex representing a dataflow operator
     *
     * @param label The operator label
     * @param inputs The operator input(s). It can be either a single channel or a list of channels.
     * @param outputs The operator output(s). It can be either a single channel, a list of channels or {@code null} if the operator has no output.
     */
    public void addOperatorNode( String label, inputs, outputs )  {
        assert label
        assert inputs
        addVertex(Type.OPERATOR, label, normalizeChannels(inputs), normalizeChannels(outputs) )
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
    void addVertex( Type type, String label, List<ChannelHandler> inbounds, List<ChannelHandler> outbounds ) {

        def vertex = createVertex( type, label)

        inbounds?.each { ChannelHandler channel ->
            inbound( vertex, channel )
        }

        outbounds?.each { ChannelHandler channel ->
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
    Vertex createVertex( Type type, String label ) {
        def result = new Vertex(type, label)
        vertices << result
        return result
    }

    private void inbound( Vertex vertex, ChannelHandler entering )  {

        def edge = findEdge(entering.instance)
        if( !edge ) {
            edges << new Edge(instance: entering.instance, to: vertex, label: entering.label)
        }
        else if( edge.to == null ) {
            edge.to = vertex
        }
        else {
            final name = getChannelName(entering)
            throw new MultipleInputChannelException(name, entering, vertex, edge.to)
        }
    }

    private void outbound( Vertex vertex, ChannelHandler leaving) {

        final edge = findEdge(leaving.instance)
        if( !edge ) {
            edges << new Edge(instance: leaving.instance, from: vertex, label: leaving.label)
        }
        else if( edge.from == null ) {
            edge.from = vertex
        }
        else {
            final name = getChannelName(leaving)
            throw new MultipleOutputChannelException(name, leaving, vertex, edge.to)
        }

    }

    static private List<ChannelHandler> normalizeInputs( InputsList inputs ) {

        inputs
                .findAll { !( it instanceof DefaultInParam)  }
                .collect { InParam p -> new ChannelHandler(instance: (DataflowChannel)p.inChannel, label: p instanceof SetInParam ? null : p.name) }

    }

    static private List<ChannelHandler> normalizeOutputs( OutputsList outputs ) {

        def result = []
        outputs.each { OutParam p ->
            if( p instanceof DefaultOutParam ) return
            p.outChannels.each {
                result << new ChannelHandler(instance: (DataflowChannel)it, label: p instanceof SetOutParam ? null : p.name)
            }
        }

        return result
    }

    static private List<ChannelHandler> normalizeChannels( entry ) {
        if( entry == null ) {
            Collections.emptyList()
        }
        else if( entry instanceof DataflowChannel ) {
            [ new ChannelHandler(instance: entry) ]
        }
        else if( entry instanceof Collection || entry instanceof Object[] ) {
            entry.collect { new ChannelHandler(instance: (DataflowChannel)it) }
        }
        else {
            throw new IllegalArgumentException("Not a valid channel type: [${entry.class.name}]")
        }
    }

    @PackageScope
    Edge findEdge( DataflowChannel channel ) {
        return edges.find { edge -> edge.instance.is(channel) }
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
    void resolveEdgeNames(Map map) {
        edges.each { Edge edge ->
            def name = resolveChannelName(map, edge.instance)
            if( name ) edge.label = name
        }
    }


    @PackageScope
    String resolveChannelName( Map map, DataflowChannel channel ) {
        def entry = map.find { k,v -> v.is channel }
        return entry ? entry.key : null
    }

    @PackageScope
    String getChannelName( ChannelHandler handler ) {
        def result = handler.label
        result ?: (session ? resolveChannelName( session.getBinding().getVariables(), handler.instance ) : null )
    }

    String render() {
        normalizeMissingVertices()
        if( session )
            resolveEdgeNames(session.getBinding().getVariables())
        else
            log.debug "Missing session object -- Cannot normalize edge names"

        def result = []
        result << "digraph graphname {"
        edges.each { edge -> result << edge.render()  }
        result << "}"
        return result.join('\n')
    }


    /**
     * Model a vertex in the DAG.
     *
     * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
     */
    @ToString
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

        /**
         * Render the DAG using the Graphviz DOT format
         * See http://www.graphviz.org/content/dot-language
         *
         * @return A string representing the DAG in the DOT notation
         */
        String render() {

            List attrs = []

            switch (type) {
                case Type.NODE:
                    attrs << "shape=point"
                    if (label) {
                        attrs << "label=\"\""
                        attrs << "xlabel=\"$label\""
                    }
                    break

                case Type.ORIGIN:
                    attrs << "shape=point"
                    attrs << "label=\"\""
                    attrs << "fixedsize=true"
                    attrs << "width=0.1"
                    if( label ) {
                        attrs << "xlabel=\"$label\""
                    }
                    break

                case Type.OPERATOR:
                    attrs << "shape=circle"
                    attrs << "label=\"\""
                    attrs << "fixedsize=true"
                    attrs << "width=0.1"
                    if( label ) {
                        attrs << "xlabel=\"$label\""
                    }
                    break

                case Type.PROCESS:
                    if( label )
                        attrs << "label=\"$label\""
                    break

                default:
                    attrs << "shape=none"
                    if( label )
                        attrs << "label=\"$label\""
            }


            return attrs ? "${getName()} [${attrs.join(',')}];" : nul
        }

    }

    /**
     * Models an edge in the DAG
     *
     * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
     */
    @PackageScope
    class Edge {

        /**
         * The {@link groovyx.gpars.dataflow.DataflowChannel} that originated this graph edge
         */
        DataflowChannel instance

        /**
         * The vertex from where the edge starts
         */
        Vertex from

        /**
         * The vertex to where the edge ends
         */
        Vertex to

        /**
         * A descriptive label
         */
        String label


        String render() {
            assert from != null || to != null

            String A = from.render()
            String B = to.render()

            def result = new StringBuilder()
            if( A ) result << A << '\n'
            if( B ) result << B << '\n'
            result << "${from.name} -> ${to.name}"
            if( label ) {
                result << " [label=\"${label}\"]"
            }
            result << ";\n"
        }

    }

    /**
     * A simple wrapper object to handle a channel and the associated label
     */
    @ToString
    static class ChannelHandler {

        /**
         * The {@link groovyx.gpars.dataflow.DataflowChannel} that originated this graph edge
         */
        DataflowChannel instance

        /**
         * The edge label
         */
        String label

    }


}