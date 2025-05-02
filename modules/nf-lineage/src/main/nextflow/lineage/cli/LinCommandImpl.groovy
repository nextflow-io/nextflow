/*
 * Copyright 2013-2025, Seqera Labs
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

import static nextflow.lineage.fs.LinPath.*

import java.nio.charset.StandardCharsets
import java.nio.file.Path

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.cli.CmdLineage
import nextflow.config.ConfigMap
import nextflow.dag.MermaidHtmlRenderer
import nextflow.lineage.LinHistoryRecord
import nextflow.lineage.LinPropertyValidator
import nextflow.lineage.LinStore
import nextflow.lineage.LinStoreFactory
import nextflow.lineage.LinUtils
import nextflow.lineage.model.FileOutput
import nextflow.lineage.model.Parameter
import nextflow.lineage.model.TaskRun
import nextflow.lineage.model.WorkflowRun
import nextflow.lineage.serde.LinEncoder
import nextflow.ui.TableBuilder
import org.eclipse.jgit.diff.DiffAlgorithm
import org.eclipse.jgit.diff.DiffFormatter
import org.eclipse.jgit.diff.RawText
import org.eclipse.jgit.diff.RawTextComparator
/**
 * Implements lineage command line operations
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class LinCommandImpl implements CmdLineage.LinCommand {

    private static final Path DEFAULT_HTML_FILE = Path.of("lineage.html")

    @Canonical
    static class Edge {
        String source
        String destination
        String label
    }

    static final private String ERR_NOT_LOADED = 'Error lineage store not loaded - Check Nextflow configuration'
    
    @Override
    void log(ConfigMap config) {
        final session = new Session(config)
        final store = LinStoreFactory.getOrCreate(session)
        if (store) {
            printHistory(store)
        } else {
            println ERR_NOT_LOADED
        }
    }

    private void printHistory(LinStore store) {
        final records = store.historyLog?.records
        if( !records ) {
            println("No workflow runs found in lineage history log")
            return
        }
        def table = new TableBuilder(cellSeparator: '\t')
            .head('TIMESTAMP')
            .head('RUN NAME')
            .head('SESSION ID')
            .head('LINEAGE ID')
        for (LinHistoryRecord record : records) {
            table.append(record.toList())
        }
        println table.toString()
    }

    @Override
    void describe(ConfigMap config, List<String> args) {
        if( !isLidUri(args[0]) )
            throw new Exception("Identifier is not a lineage URL")
        final store = LinStoreFactory.getOrCreate(new Session(config))
        if ( !store ) {
            println ERR_NOT_LOADED
            return
        }
        try {
            def entry = LinUtils.getMetadataObject(store, new URI(args[0]))
            if( !entry ) {
                println "No entry found for ${args[0]}"
                return
            }
            println LinUtils.encodeSearchOutputs(entry, true)
        } catch (Throwable e) {
            println "Error loading ${args[0]} - ${e.message}"
        }
    }

    @Override
    void render(ConfigMap config, List<String> args) {
        final store = LinStoreFactory.getOrCreate(new Session(config))
        if( !store ) {
            println ERR_NOT_LOADED
            return
        }
        try {
            final renderFile = args.size() > 1 ? Path.of(args[1]) : DEFAULT_HTML_FILE
            renderLineage(store, args[0], renderFile)
            println("Linage graph for ${args[0]} rendered in $renderFile")
        } catch (Throwable e) {
            println("ERROR: rendering lineage graph - ${e.message}")
        }
    }

    private void renderLineage(LinStore store, String dataLid, Path file) {
        def lines = [] as List<String>
        lines << "flowchart BT".toString()
        final nodesToRender = new LinkedList<String>()
        nodesToRender.add(dataLid)
        final edgesToRender = new LinkedList<Edge>()
        while (!nodesToRender.isEmpty()) {
            final node = nodesToRender.removeFirst()
            processNode(lines, node, nodesToRender, edgesToRender, store)
        }
        lines << ""
        edgesToRender.each { lines << "    ${it.source} -->${it.destination}".toString() }
        lines << ""
        lines.join('\n')
        final template = MermaidHtmlRenderer.readTemplate()
        file.text = template.replace('REPLACE_WITH_NETWORK_DATA', lines.join('\n'))
    }

    private String safeId( String rawId){
        return rawId.replaceAll(/[^a-zA-Z0-9_.:\/\-]/, '_')
    }

    private void processNode(List<String> lines, String nodeToRender, LinkedList<String> nodes, LinkedList<Edge> edges, LinStore store) {
        if (!isLidUri(nodeToRender))
            throw new Exception("Identifier is not a lineage URL")
        final key = nodeToRender.substring(LID_PROT.size())
        final lidObject = store.load(key)
        switch (lidObject.getClass()) {
            case FileOutput:
                processDataOutput(lidObject as FileOutput, lines, nodeToRender, nodes, edges)
                break;

            case WorkflowRun:
                processWorkflowRun(lidObject as WorkflowRun, lines, nodeToRender, edges)
                break

            case TaskRun:
                processTaskRun(lidObject as TaskRun, lines, nodeToRender, nodes, edges)
                break

            default:
                throw new Exception("Unrecognized type reference ${lidObject.getClass().getSimpleName()}")
        }
    }

    private void processTaskRun(TaskRun taskRun, List<String> lines, String nodeToRender, LinkedList<String> nodes, LinkedList<Edge> edges) {
        lines << "    ${nodeToRender}@{shape: process, label: \"${taskRun.name} [$nodeToRender]\"}".toString()
        final parameters = taskRun.input
        for (Parameter source : parameters) {
            if (source.type.equals("path")) {
                manageFileInParam(lines, nodeToRender, nodes, edges, source.value)
            } else {
                final label = convertToLabel(source.value.toString())
                final id = safeId(source.value.toString())
                lines << "    ${id}@{shape: document, label: \"${label}\"}".toString();
                edges.add(new Edge(id, nodeToRender))
            }
        }
    }

    private void processWorkflowRun(WorkflowRun wfRun, List<String> lines, String nodeToRender, LinkedList<Edge> edges) {
        lines << """    ${nodeToRender}@{shape: processes, label: \"${wfRun.name} [${nodeToRender}]\"}""".toString()
        final parameters = wfRun.params
        parameters.each {
            final label = convertToLabel(it.value.toString())
            final id = safeId(it.value.toString())
            lines << "    ${id}@{shape: document, label: \"${label}\"}".toString();
            edges.add(new Edge(id, nodeToRender))
        }
    }

    private void processDataOutput(FileOutput lidObject, List<String> lines, String nodeToRender, LinkedList<String> nodes, LinkedList<Edge> edges){
        lines << "    ${nodeToRender}@{shape: document, label: \"${nodeToRender}\"}".toString();
        final source = lidObject.source
        if(! source )
            return
        if (isLidUri(source)) {
            nodes.add(source)
            edges.add(new Edge(source, nodeToRender))
        } else {
            final label = convertToLabel(source)
            final id = safeId(source)
            lines << "    ${id}@{shape: document, label: \"${label}\"}".toString();
            edges.add(new Edge(id, nodeToRender))
        }
    }

    private String convertToLabel(String label){
        return label.replace('http', 'h\u200Ettp')
    }

    private void manageFileInParam(List<String> lines, String nodeToRender, LinkedList<String> nodes, LinkedList<Edge> edges, value){
        if (value instanceof Collection) {
            value.each { manageFileInParam(lines, nodeToRender, nodes, edges, it) }
            return
        }
        if (value instanceof CharSequence) {
            final source = value.toString()
            if (isLidUri(source)) {
                nodes.add(source)
                edges.add(new Edge(source, nodeToRender))
                return
            }
        }
        if (value instanceof Map ) {
            if (value.path) {
                final path = value.path.toString()
                if (isLidUri(path)) {
                    nodes.add(path)
                    edges.add(new Edge(path, nodeToRender))
                    return
                } else {
                    final label = convertToLabel(path)
                    final id = safeId(path)
                    lines << "    ${id}@{shape: document, label: \"${label}\"}".toString();
                    edges.add(new Edge(id, nodeToRender))
                    return
                }
            }
        }
        final label = convertToLabel(value.toString())
        final id = safeId(value.toString())
        lines << "    ${id}@{shape: document, label: \"${label}\"}".toString();
        edges.add(new Edge(id, nodeToRender))
    }

    @Override
    void diff(ConfigMap config, List<String> args) {
        if (!isLidUri(args[0]) || !isLidUri(args[1]))
            throw new Exception("Identifier is not a lineage URL")

        final store = LinStoreFactory.getOrCreate(new Session(config))
        if (!store) {
            println ERR_NOT_LOADED
            return
        }
        try {
            final key1 = args[0].substring(LID_PROT.size())
            final entry1 = store.load(key1)
            if (!entry1) {
                println "No entry found for ${args[0]}."
                return
            }
            final key2 = args[1].substring(LID_PROT.size())
            final entry2 = store.load(key2)
            if (!entry2) {
                println "No entry found for ${args[1]}."
                return
            }
            final encoder = new LinEncoder().withPrettyPrint(true)
            generateDiff(encoder.encode(entry1), key1, encoder.encode(entry2), key2)
        } catch (Throwable e) {
            println "Error generating diff between ${args[0]}: $e.message"
        }
    }

    private static void generateDiff(String entry1, String key1, String entry2, String key2) {
        // Convert strings to JGit RawText format
        final text1 = new RawText(entry1.getBytes(StandardCharsets.UTF_8))
        final text2 = new RawText(entry2.getBytes(StandardCharsets.UTF_8))

        // Set up the diff algorithm (Git-style diff)
        final diffAlgorithm = DiffAlgorithm.getAlgorithm(DiffAlgorithm.SupportedAlgorithm.MYERS)
        final diffComparator = RawTextComparator.DEFAULT

        // Compute the differences
        final editList = diffAlgorithm.diff(diffComparator, text1, text2)

        final output = new StringBuilder()
        // Add header
        output.append("diff --git ${key1} ${key2}\n")
        output.append("--- ${key1}\n")
        output.append("+++ ${key2}\n")

        // Use DiffFormatter to display results in Git-style format
        final outputStream = new ByteArrayOutputStream()
        final diffFormatter = new DiffFormatter(outputStream)
        diffFormatter.setOldPrefix(key1)
        diffFormatter.setNewPrefix(key2)
        diffFormatter.format(editList, text1, text2)
        output.append(outputStream.toString(StandardCharsets.UTF_8))

        println output.toString()
    }

    @Override
    void find(ConfigMap config, List<String> args) {
        final store = LinStoreFactory.getOrCreate(new Session(config))
        if (!store) {
            println ERR_NOT_LOADED
            return
        }
        try {
            final params = parseFindArgs(args)
            new LinPropertyValidator().validateQueryParams(params)
            println LinUtils.encodeSearchOutputs( store.search(params).keySet().collect { asUriString(it) }, true )
        } catch (Throwable e){
            println "Error searching for ${args[0]}. ${e.message}"
        }
    }

    private Map<String, List<String>> parseFindArgs(List<String> args){
        Map<String, List<String>> params = [:].withDefault { [] }

        args.collectEntries { pair ->
            def idx = pair.indexOf('=')
            if( idx < 0 )
                throw new IllegalArgumentException("Parameter $pair doesn't contain '=' separator")
            final key = URLDecoder.decode(pair[0..<idx], 'UTF-8')
            final value = URLDecoder.decode(pair[(idx + 1)..<pair.length()], 'UTF-8')
            params[key] << value
        }
        return params
    }
}
