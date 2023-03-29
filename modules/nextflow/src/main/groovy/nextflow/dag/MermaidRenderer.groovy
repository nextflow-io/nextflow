/*
 * Copyright 2013-2023, Seqera Labs
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
    void renderAbstractGraph(DAG dag, Path file) {
        def lines = []
        lines << "flowchart TD"

        dag.vertices.each { vertex ->
            lines << "    ${renderVertex( vertex )}"
        }

        dag.edges.each { edge ->
            lines << "    ${renderEdge( edge )}"
        }

        lines << ""

        file.text = lines.join('\n')
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

    @Override
    void renderConcreteGraph(ConcreteDAG graph, Path file) {
        def renderedOutputs = [] as Set<Path>
        def numInputs = 0
        def numOutputs = 0

        def lines = []
        lines << "flowchart TD"

        // render tasks and task inputs
        graph.nodes.values().each { task ->
            // render task node
            lines << "    ${task.getSlug()}[\"${task.label}\"]"

            task.inputs.each { input ->
                // render task input from predecessor
                if( input.predecessor != null ) {
                    final pred = graph.nodes[input.predecessor]
                    lines << "    ${pred.getSlug()} -->|${input.name}| ${task.getSlug()}"
                    renderedOutputs << input.path
                }

                // render task input from source node
                else {
                    numInputs += 1
                    lines << "    i${numInputs}(( )) -->|${input.name}| ${task.getSlug()}"
                }
            }
        }

        // render task outputs with sink nodes
        graph.nodes.values().each { task ->
            task.outputs.each { output ->
                if( output.path !in renderedOutputs ) {
                    numOutputs += 1
                    lines << "    ${task.getSlug()} -->|${output.name}| o${numOutputs}(( ))"
                }
            }
        }

        lines << ""

        file.text = lines.join('\n')
    }
}
