package nextflow.dag

import groovy.transform.CompileStatic
import java.nio.file.Path

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Mike Smoot <mes@aescon.com>
 */
class DagreD3Renderer implements DagRenderer {

    /**
     * Render the DAG in HTML using Cytoscape.js
     * to the specified file.
     * See http://js.cytoscape.org for more info.
     */
    @Override
    void renderDocument(DAG dag, Path file) {
        String tmplPage = readTemplate()
        String network = renderNetwork(dag)
        file.text = tmplPage.replaceAll(~/\/\* REPLACE_WITH_NETWORK_DATA \*\//, network)
    }

    private String readTemplate() {
        StringWriter writer = new StringWriter();
        def res = DagreD3Renderer.class.getResourceAsStream('dagre.d3.dag.template.html')
        int ch
        while( (ch=res.read()) != -1 ) {
            writer.append(ch as char);
        }
        writer.toString();
    }

    static String renderNetwork(DAG dag) {
        def result = []

        dag.vertices.each { vertex -> result << renderVertex( vertex ) }
        dag.edges.each { edge -> result << renderEdge( edge ) }

        return result.join('\n')
    }

    private static String renderVertex(vertex) {
        String dat = "g.setNode('${vertex.getName()}'"
        if (vertex.label) {
            return dat + ", { label: '${vertex.label}' });"
        } else {
            return dat + ", { label: '' });"
        }
    }


    private static String renderEdge(edge) {
        assert edge.from != null && edge.to != null
        String dat = "g.setEdge('${edge.from.name}', '${edge.to.name}'"
        if ( edge.label ) {
            return dat + ", { label: '${edge.label}' });"
        }
        else {
            return dat + ", { label: '' });"
        }
    }

}
