/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.lineage.cli

import java.nio.file.Path

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.dag.MermaidHtmlRenderer
import nextflow.lineage.LinStore
import nextflow.lineage.model.v1beta1.FileOutput
import nextflow.lineage.model.v1beta1.TaskRun
import nextflow.lineage.model.v1beta1.WorkflowRun

import static nextflow.lineage.fs.LinPath.LID_PROT
import static nextflow.lineage.fs.LinPath.isLidUri
/**
 * Renderer for the lineage graph.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class LinDagRenderer {

    private final LinStore store

    private Queue<String> queue

    private List<Node> nodes

    private List<Edge> edges

    LinDagRenderer(LinStore store) {
        this.store = store
    }

    private void enqueueLid(String lid) {
        queue.add(lid)
    }

    private void addNode(String id, String label, NodeType type) {
        nodes.add(new Node(id, label, type))
    }

    private void addEdge(String source, String target) {
        edges.add(new Edge(source, target))
    }

    void render(String lid, Path file) {
        // visit nodes
        this.queue = new LinkedList<String>()
        this.nodes = new LinkedList<Node>()
        this.edges = new LinkedList<Edge>()

        enqueueLid(lid)
        while( !queue.isEmpty() )
            visitLid(queue.remove())

        // render Mermaid diagram
        final lines = new LinkedList<String>()
        lines << "flowchart TB"

        for( final node : nodes ) {
            if( node.type == NodeType.FILE )
                lines << "    ${node.id}[\"${node.label}\"]".toString()
            if( node.type == NodeType.TASK )
                lines << "    ${node.id}([\"${node.label}\"])".toString()
        }

        for( final edge : edges )
            lines << "    ${edge.source} --> ${edge.target}".toString()

        final template = MermaidHtmlRenderer.readTemplate()
        file.text = template.replace('REPLACE_WITH_NETWORK_DATA', lines.join('\n'))
    }

    private void visitLid(String lid) {
        if( !isLidUri(lid) )
            throw new Exception("Identifier is not a lineage URL: ${lid}")
        final record = store.load(rawLid(lid))
        if( !record ) {
            log.warn "Lineage record references an LID that does not exist: ${lid}"
            return
        }
        if( record instanceof FileOutput )
            visitFileOutput(lid, record)
        else if( record instanceof TaskRun )
            visitTaskRun(lid, record)
        else if( record instanceof WorkflowRun )
            visitWorkflowRun(lid, record)
        else
            throw new Exception("Cannot render lineage for type ${record.getClass().getSimpleName()} -- must be a FileOutput, TaskRun, or WorkflowRun")
    }

    private void visitFileOutput(String lid, FileOutput fileOutput) {
        addNode(lid, lid, NodeType.FILE)
        final source = fileOutput.source
        if( !source )
            return
        if( isLidUri(source) ) {
            enqueueLid(source)
            addEdge(source, lid)
        }
        else {
            final id = safeId(source)
            final label = safeLabel(source)
            addNode(id, label, NodeType.FILE)
            addEdge(id, lid)
        }
    }

    private void visitTaskRun(String lid, TaskRun taskRun) {
        addNode(lid, "${taskRun.name} [${lid}]", NodeType.TASK)
        for( final param : taskRun.input ) {
            visitParameter(lid, param.value)
        }
    }

    private void visitWorkflowRun(String lid, WorkflowRun workflowRun) {
        addNode(lid, "${workflowRun.name} [${lid}]", NodeType.TASK)
        for( final param : workflowRun.params ) {
            visitParameter0(lid, param.value.toString())
        }
    }

    private void visitParameter(String lid, Object value) {
        if( value instanceof Collection ) {
            for( final el : value )
                visitParameter(lid, el)
        }
        else if( value instanceof CharSequence ) {
            final source = value.toString()
            if( isLidUri(source) ) {
                enqueueLid(source)
                addEdge(source, lid)
            }
            else {
                visitParameter0(lid, source)
            }
        }
        else if( value instanceof Map && value.path ) {
            final path = value.path.toString()
            if( isLidUri(path) ) {
                enqueueLid(path)
                addEdge(path, lid)
            }
            else {
                visitParameter0(lid, path)
            }
        }
        else {
            visitParameter0(lid, value.toString())
        }
    }

    private void visitParameter0(String lid, String value) {
        final id = safeId(value)
        final label = safeLabel(value)
        addNode(id, label, NodeType.FILE)
        addEdge(id, lid)
    }

    private static String rawLid(String lid) {
        return lid.substring(LID_PROT.size())
    }

    private static String safeId(String rawId) {
        return rawId.replaceAll(/[^a-zA-Z0-9_.:\/\-]/, '_')
    }

    private static String safeLabel(String label) {
        return label.replace('http', 'h\u200Ettp')
    }

    @Canonical
    private static class Node {
        String id
        String label
        NodeType type
    }

    @Canonical
    private static class Edge {
        String source
        String target
    }

    private static enum NodeType {
        FILE,
        TASK,
    }

}
