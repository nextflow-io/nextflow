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
import java.nio.file.Path

/**
 * Render the DAG in Cytoscape.js compatible
 * JSON to the specified file.
 * See http://js.cytoscape.org for more info.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Mike Smoot <mes@aescon.com>
 */
class CytoscapeJsRenderer implements DagRenderer {

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
