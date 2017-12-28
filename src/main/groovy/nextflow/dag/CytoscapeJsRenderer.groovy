/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
