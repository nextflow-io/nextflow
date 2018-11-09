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

import groovy.transform.PackageScope

/**
 * Render the DAG using the Graphviz DOT format
 * to the specified file.
 * See http://www.graphviz.org/content/dot-language
 * for more info.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Mike Smoot <mes@aescon.com>
 */
class DotRenderer implements DagRenderer {

    private String name

    /**
     * Create a render instance
     *
     * @param name The graph name used in the DOT format
     */
    DotRenderer( String name ) {
        this.name = normalise(name)
    }

    @PackageScope
    static String normalise(String str) { str.replaceAll(/[^0-9_A-Za-z]/,'') }

    @Override
    void renderDocument(DAG dag, Path file) {
        file.text = renderNetwork(dag)
    }

    String renderNetwork(DAG dag) {
        def result = []
        result << "digraph \"$name\" {"
        dag.edges.each { edge -> result << renderEdge( edge ) }
        result << "}\n"
        return result.join('\n')
    }

    private static String renderVertex(DAG.Vertex vertex) {

        List attrs = []

        switch (vertex.type) {
            case DAG.Type.NODE:
                attrs << "shape=point"
                if (vertex.label) {
                    attrs << "label=\"\""
                    attrs << "xlabel=\"$vertex.label\""
                }
                break

            case DAG.Type.ORIGIN:
                attrs << "shape=point"
                attrs << "label=\"\""
                attrs << "fixedsize=true"
                attrs << "width=0.1"
                if( vertex.label ) {
                    attrs << "xlabel=\"$vertex.label\""
                }
                break

            case DAG.Type.OPERATOR:
                attrs << "shape=circle"
                attrs << "label=\"\""
                attrs << "fixedsize=true"
                attrs << "width=0.1"
                if( vertex.label ) {
                    attrs << "xlabel=\"$vertex.label\""
                }
                break

            case DAG.Type.PROCESS:
                if( vertex.label )
                    attrs << "label=\"$vertex.label\""
                break

            default:
                attrs << "shape=none"
                if( vertex.label )
                    attrs << "label=\"$vertex.label\""
        }


        return attrs ? "${vertex.getName()} [${attrs.join(',')}];" : null
    }

    private static String renderEdge(DAG.Edge edge) {
        assert edge.from != null && edge.to != null

        String A = renderVertex( edge.from )
        String B = renderVertex( edge.to )

        def result = new StringBuilder()
        if( A ) result << A << '\n'
        if( B ) result << B << '\n'
        result << "${edge.from.name} -> ${edge.to.name}"
        if( edge.label ) {
            result << " [label=\"${edge.label}\"]"
        }
        result << ";\n"
    }

}
