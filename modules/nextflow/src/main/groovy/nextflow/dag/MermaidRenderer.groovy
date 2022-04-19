/*
 * Copyright 2020-2022, Seqera Labs
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
 * Render the DAG using the Mermaid format to the specified file.
 * See https://mermaid-js.github.io/mermaid/#/ for more info.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class MermaidRenderer implements DagRenderer {

    @Override
    void renderDocument(DAG dag, Path file) {
        file.text = renderNetwork(dag)
    }

    String renderNetwork(DAG dag) {
        def lines = []
        lines << "flowchart TD"

        dag.vertices.each { vertex ->
            lines << "    ${renderVertex( vertex )}"
        }

        dag.edges.each { edge ->
            lines << "    ${renderEdge( edge )}"
        }

        lines << ""

        return lines.join('\n')
    }

    private String renderVertex(DAG.Vertex vertex) {
        final id = vertex.getName()

        switch (vertex.type) {
            case DAG.Type.NODE:
                return "${id}((${vertex.label ?: ' '}))"

            case DAG.Type.ORIGIN:
                return "${id}((${vertex.label ?: ' '}))"

            case DAG.Type.OPERATOR:
                return "${id}([${vertex.label ?: ' '}])"

            case DAG.Type.PROCESS:
                return "${id}[${vertex.label ?: ' '}]"

            default:
                return "${id}[${vertex.label ?: ' '}]"
        }
    }

    private String renderEdge(DAG.Edge edge) {
        assert edge.from != null && edge.to != null

        String label = edge.label ? "|${edge.label}|" : ""

        return "${edge.from.name} -->${label} ${edge.to.name}"
    }
}
