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

import groovy.transform.CompileStatic
import nextflow.Global
import nextflow.Session
import nextflow.script.WorkflowMetadata

/**
 * Render the DAG using the Mermaid format to the specified file.
 * See https://mermaid-js.github.io/mermaid/#/ for more info.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class MermaidRenderer implements DagRenderer {

    private WorkflowMetadata metadata = ((Session)Global.session).workflowMetadata

    @Override
    void renderWorkflowGraph(DAG dag, Path file) {
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
    void renderTaskGraph(TaskDAG dag, Path file) {
        def lines = []
        lines << "flowchart TD"

        // render workflow inputs
        final inputs = [:] as Map<Path,String>

        lines << "    subgraph \" \""

        dag.vertices.each { task, vertex ->
            vertex.inputs.each { name, path ->
                if( dag.getProducerVertex(path) || path in inputs )
                    return

                inputs[path] = "in${inputs.size()}".toString()
                lines << "    ${inputs[path]}[\"${normalizePath(path)}\"]"
            }
        }

        lines << "    end"

        // render tasks
        final taskOutputs = [] as Set<Path>

        dag.vertices.each { task, vertex ->
            lines << "    ${vertex.getSlug()}([\"${vertex.label}\"])"

            vertex.inputs.each { name, path ->
                // render task input from predecessor task
                final pred = dag.getProducerVertex(path)
                if( pred != null ) {
                    lines << "    ${pred.getSlug()} -->|${path.name}| ${vertex.getSlug()}"
                    taskOutputs << path
                }

                // render task input from source node
                else {
                    lines << "    ${inputs[path]} --> ${vertex.getSlug()}"
                }
            }
        }

        // render task outputs
        final outputs = [:] as Map<Path,String>

        dag.vertices.each { task, vertex ->
            vertex.outputs.each { path ->
                if( path in taskOutputs || path in outputs )
                    return

                outputs[path] = "out${outputs.size()}".toString()
                lines << "    ${vertex.getSlug()} --> ${outputs[path]}"
            }
        }

        // render workflow outputs
        lines << "    subgraph \" \""

        outputs.each { path, label ->
            lines << "    ${label}[${path.name}]"
        }

        lines << "    end"

        lines << ""

        file.text = lines.join('\n')
    }

    String normalizePath(Path path) {
        path.toUriString()
            .replace(metadata.workDir.toUriString(), '/work')
            .replace(metadata.projectDir.toUriString(), '')
            .replace(metadata.launchDir.toUriString(), '')
    }
}
