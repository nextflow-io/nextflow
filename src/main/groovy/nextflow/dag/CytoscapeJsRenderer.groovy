package nextflow.dag

import groovy.transform.CompileStatic
import java.nio.file.Path

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Mike Smoot <mes@aescon.com>
 */
class CytoscapeJsRenderer implements DagRenderer {

    /**
     * Render the DAG in Cytoscape.js compatibile
     * JSON to the specified file.
     * See http://js.cytoscape.org for more info.
     */
    @Override
    void renderDocument(DAG dag, Path file) {
        file.text = renderNetwork(dag)
    }

    static String renderNetwork(DAG dag) {
        def result = []
        result << "elements: {"

        result << "nodes: ["
        dag.vertices.each { vertex -> result << renderVertex( vertex ) }
        result << "],"

        result << "edges: ["
        dag.edges.each { edge -> result << renderEdge( edge ) }
        result << "],"

        result << "},"

        return result.join('\n')
    }

    private static String renderVertex(vertex) {
        String pre = "{ data: { id: '${vertex.getName()}'"
        String post = "}, classes: '${vertex.type.name()}' },"
        if (vertex.label) {
            return pre + ", label: '${vertex.label}'" + post
        } else {
            return pre + post
        }
    }


    private static String renderEdge(edge) {
        assert edge.from != null && edge.to != null
        String dat = "{ data: { source: '${edge.from.name}', target: '${edge.to.name}'"
        if ( edge.label ) {
            return dat + ", label: '${edge.label}' } },"
        }
        else {
            return dat + "} },"
        }
    }

}
