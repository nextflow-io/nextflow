/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

/**
 * Render the DAG in HTML using Cytoscape.js
 * to the specified file.
 * See http://js.cytoscape.org for more info.
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Mike Smoot <mes@aescon.com>
 */
class CytoscapeHtmlRenderer implements DagRenderer {

    @Override
    void renderDocument(DAG dag, Path file) {
        String tmplPage = readTemplate()
        String network = CytoscapeJsRenderer.renderNetwork(dag)
        file.text = tmplPage.replaceAll(~/\/\* REPLACE_WITH_NETWORK_DATA \*\//, network)
    }

    private String readTemplate() {
        StringWriter writer = new StringWriter();
        def res = CytoscapeHtmlRenderer.class.getResourceAsStream('cytoscape.js.dag.template.html')
        int ch
        while( (ch=res.read()) != -1 ) {
            writer.append(ch as char);
        }
        writer.toString();
    }
}
