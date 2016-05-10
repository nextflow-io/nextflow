package nextflow.dag

import groovy.transform.CompileStatic
import java.nio.file.Path

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Mike Smoot <mes@aescon.com>
 */
class CytoscapeHtmlRenderer implements DagRenderer {

    /**
     * Render the DAG in HTML using Cytoscape.js
     * to the specified file.
     * See http://js.cytoscape.org for more info.
     */
    @Override
    void renderDocument(DAG dag, Path file) {
        String tmplPage = readTemplate()
        String network = CytoscapeJsRenderer.renderNetwork(dag)
        file.text = tmplPage.replaceAll(~/\/\* REPLACE_WITH_NETWORK_DATA \*\//, network)
    }

    private String readTemplate() {
        StringWriter writer = new StringWriter();
        def res = CytoscapeJsRenderer.class.getResourceAsStream('cytoscape.js.dag.template.html')
        int ch
        while( (ch=res.read()) != -1 ) {
            writer.append(ch as char);
        }
        writer.toString();
    }
}
