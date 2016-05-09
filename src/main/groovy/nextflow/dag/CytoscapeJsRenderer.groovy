package nextflow.dag

import groovy.transform.PackageScope
import groovy.transform.ToString
import nextflow.Session


/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Mike Smoot <mes@aescon.com>
 */
class CytoscapeJsRenderer {

    /**
     * Render the DAG in HTML using Cytoscape.js
     * See http://js.cytoscape.org
     *
     * @return A string representing the DAG in
     *         HTML rendered using Cytoscape.js.
     */
    static String render(DAG dag) {
        String tmplPage = readTemplate()
        String network = renderNetwork(dag)
        return tmplPage.replaceAll(~/\/\* REPLACE_WITH_NETWORK_DATA \*\//, network)
    }

    private static String readTemplate() {
        StringWriter writer = new StringWriter();
        def res = CytoscapeJsRenderer.class.getResourceAsStream('cytoscape.js.dag.template.html')
        int ch
        while( (ch=res.read()) != -1 ) {
            writer.append(ch as char);
        }
        writer.toString();
    }

    private static String renderNetwork(DAG dag) {
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
