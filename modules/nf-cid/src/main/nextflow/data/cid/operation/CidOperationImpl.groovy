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
 *
 */

package nextflow.data.cid.operation

import org.eclipse.jgit.diff.DiffAlgorithm
import org.eclipse.jgit.diff.DiffFormatter
import org.eclipse.jgit.diff.RawText
import org.eclipse.jgit.diff.RawTextComparator

import java.nio.charset.StandardCharsets

import static nextflow.data.cid.fs.CidPath.*

import java.nio.file.Path

import groovy.json.JsonSlurper
import groovy.transform.Canonical
import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.cli.CmdCid
import nextflow.config.ConfigMap
import nextflow.dag.MermaidHtmlRenderer
import nextflow.data.cid.CidHistoryRecord
import nextflow.data.cid.CidStore
import nextflow.data.cid.CidStoreFactory
import nextflow.data.cid.model.DataType
import nextflow.ui.TableBuilder
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CidOperationImpl implements CmdCid.CidOperation {

    @Canonical
    static class Edge {
        String source
        String destination
        String label
    }

    @Override
    void log(ConfigMap config) {
        final session = new Session(config)
        final store = CidStoreFactory.getOrCreate(session)
        if (store) {
            printHistory(store)
        } else {
            println "Error CID store not loaded. Check Nextflow configuration."
        }
    }

    private void printHistory(CidStore store) {
        final records = store.historyLog?.records
        if( records ) {
            def table = new TableBuilder(cellSeparator: '\t')
                .head('TIMESTAMP')
                .head('RUN NAME')
                .head('SESSION ID')
                .head('RUN CID')
                .head('RESULT CID')
            for( CidHistoryRecord record: records ){
                table.append(record.toList())
            }
            println table.toString()
        } else {
            println("No workflow runs CIDs found.")
        }
    }

    @Override
    void show(ConfigMap config, List<String> args) {
        if (!args[0].startsWith(CID_PROT))
            throw new Exception("Identifier is not a CID URL")
        final key = args[0].substring(CID_PROT.size())
        final store = CidStoreFactory.getOrCreate(new Session(config))
        if (store) {
            try {
                final entry = store.load(key)
                if( entry )
                    println entry.toString()
                else
                    println "No entry found for ${args[0]}."
            } catch (Throwable e) {
                println "Error loading ${args[0]}."
            }
        } else {
            println "Error CID store not loaded. Check Nextflow configuration."
        }
    }

    @Override
    void lineage(ConfigMap config, List<String> args) {
        try {
            final store = CidStoreFactory.getOrCreate(new Session(config))
            final template = MermaidHtmlRenderer.readTemplate()
            final network = getLineage(store, args[0])
            Path file = Path.of(args[1])
            file.text = template.replace('REPLACE_WITH_NETWORK_DATA', network)
            println("Linage graph for ${args[0]} rendered in ${args[1]}")
        } catch (Throwable e) {
            println("ERROR: rendering lineage graph. ${e.message}")
        }
    }

    private String getLineage(CidStore store, String dataCid) {
        def lines = [] as List<String>
        lines << "flowchart BT".toString()

        final nodesToRender = new LinkedList<String>()
        nodesToRender.add(dataCid)
        final edgesToRender = new LinkedList<Edge>()
        while (!nodesToRender.isEmpty()) {
            final node = nodesToRender.removeFirst()
            processNode(lines, node, nodesToRender, edgesToRender, store)
        }
        lines << ""
        edgesToRender.each { lines << "    ${it.source} -->${it.destination}".toString() }
        lines << ""
        return lines.join('\n')
    }

    private void processNode(List<String> lines, String nodeToRender, LinkedList<String> nodes, LinkedList<Edge> edges, CidStore store) {
        if (!nodeToRender.startsWith(CID_PROT))
            throw new Exception("Identifier is not a CID URL")
        final slurper = new JsonSlurper()
        final key = nodeToRender.substring(CID_PROT.size())
        final cidObject = slurper.parse(store.load(key).toString().toCharArray()) as Map
        switch (DataType.valueOf(cidObject.type as String)) {
            case DataType.TaskOutput:
            case DataType.WorkflowOutput:
                lines << "    ${nodeToRender}@{shape: document, label: \"${nodeToRender}\"}".toString();
                final source = cidObject.source as String
                if (source) {
                    if (source.startsWith(CID_PROT)) {
                        nodes.add(source)
                        edges.add(new Edge(source, nodeToRender))
                    } else {
                        final label = convertToLabel(source)
                        lines << "    ${source}@{shape: document, label: \"${label}\"}".toString();
                        edges.add(new Edge(source, nodeToRender))
                    }
                }

                break;
            case DataType.WorkflowRun:
                lines << "${nodeToRender}@{shape: processes, label: \"${cidObject.runName}\"}".toString()
                final parameters = cidObject.params as List<nextflow.data.cid.model.Parameter>
                parameters.each {
                    final label = convertToLabel(it.value.toString())
                    lines << "    ${it.value.toString()}@{shape: document, label: \"${label}\"}".toString();
                    edges.add(new Edge(it.value.toString(), nodeToRender))
                }
                break;
            case DataType.TaskRun:
                lines << "    ${nodeToRender}@{shape: process, label: \"${cidObject.name}\"}".toString()
                final parameters = cidObject.inputs as List<nextflow.data.cid.model.Parameter>
                for (nextflow.data.cid.model.Parameter source: parameters){
                    if (source.type.equals(nextflow.script.params.FileInParam.simpleName)) {
                        manageFileInParam(lines, nodeToRender, nodes, edges, source.value)
                    } else {
                        final label = convertToLabel(source.value.toString())
                        lines << "    ${source.value.toString()}@{shape: document, label: \"${label}\"}".toString();
                        edges.add(new Edge(source.value.toString(), nodeToRender))
                    }
                }
                break;
            default:
                throw new Exception("Unrecognized type reference ${cidObject.type}")
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
            if (source.startsWith(CID_PROT)) {
                nodes.add(source)
                edges.add(new Edge(source, nodeToRender))
                return
            }
        }
        if (value instanceof Map) {
            if (value.path) {
                final path = value.path.toString()
                if (path.startsWith(CID_PROT)) {
                    nodes.add(path)
                    edges.add(new Edge(path, nodeToRender))
                    return
                } else {
                    final label = convertToLabel(path)
                    lines << "    ${path}@{shape: document, label: \"${label}\"}".toString();
                    edges.add(new Edge(path, nodeToRender))
                    return
                }
            }
        }
        final label = convertToLabel(value.toString())
        lines << "    ${value.toString()}@{shape: document, label: \"${label}\"}".toString();
        edges.add(new Edge(value.toString(), nodeToRender))
    }

    @Override
    void diff(ConfigMap config, List<String> args) {
        if (!args[0].startsWith(CID_PROT) || !args[1].startsWith(CID_PROT))
            throw new Exception("Identifier is not a CID URL")

        final store = CidStoreFactory.getOrCreate(new Session(config))
        if (store) {
            try {
                final key1 = args[0].substring(CID_PROT.size())
                final entry1 = store.load(key1) as String
                if( !entry1 ){
                    println "No entry found for ${args[0]}."
                    return
                }
                final key2 = args[1].substring(CID_PROT.size())
                final entry2 = store.load(key2) as String
                if( !entry2 ) {
                    println "No entry found for ${args[1]}."
                    return
                }
                generateDiff(entry1, key1, entry2, key2)
            } catch (Throwable e) {
                println "Error generating diff between ${args[0]}: $e.message"
            }
        } else {
            println "Error CID store not loaded. Check Nextflow configuration."
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


}
