/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.cli

import com.beust.jcommander.Parameter
import groovy.json.JsonSlurper
import groovy.transform.Canonical
import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.config.ConfigBuilder
import nextflow.dag.MermaidHtmlRenderer
import nextflow.data.cid.CidHistoryFile
import nextflow.data.cid.CidStore
import nextflow.data.cid.CidStoreFactory
import nextflow.data.cid.model.DataType
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.ui.TableBuilder

import java.nio.file.Path
import java.nio.file.Paths

import static nextflow.data.cid.fs.CidPath.CID_PROT
import static nextflow.data.cid.fs.CidPath.METADATA_FILE

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CmdCid extends CmdBase {

    private static final String NAME = 'cid'

    interface SubCmd {
        String getName()
        void apply(List<String> args)
        void usage()
    }

    private List<SubCmd> commands = new ArrayList<>()

    CmdCid() {
        commands << new CmdLog()
        commands << new CmdShow()
        commands << new CmdLineage()
    }

    @Parameter(hidden = true)
    List<String> args

    @Override
    String getName() {
        return NAME
    }

    @Override
    void run() {
        if( !args ) {
            return
        }
        // setup the plugins system and load the secrets provider
        Plugins.init()

        getCmd(args).apply(args.drop(1))
    }

    protected SubCmd getCmd(List<String> args) {

        def cmd = commands.find { it.name == args[0] }
        if( cmd ) {
            return cmd
        }

        def matches = commands.collect{ it.name }.closest(args[0])
        def msg = "Unknown cloud sub-command: ${args[0]}"
        if( matches )
            msg += " -- Did you mean one of these?\n" + matches.collect { "  $it"}.join('\n')
        throw new AbortOperationException(msg)
    }

    class CmdLog implements SubCmd {

        @Override
        String getName() {
            return 'log'
        }

        @Override
        void apply(List<String> args) {
            if (args.size() != 0) {
                println("ERROR: Incorrect number of parameters")
                usage()
                return
            }
            final config = new ConfigBuilder()
                    .setOptions(getLauncher().getOptions())
                    .setBaseDir(Paths.get('.'))
                    .build()
            final session = new Session(config)
            final store = CidStoreFactory.getOrCreate(session)
            if (store) {
                printHistory(store)
            } else {
                println "Error CID store not loaded. Check Nextflow configuration."
            }
        }

        private void printHistory(CidStore store) {
            final historyFile = store.getHistoryFile()
            if (historyFile.exists()) {
                def table = new TableBuilder(cellSeparator: '\t')
                    .head('TIMESTAMP')
                    .head('RUN NAME')
                    .head('SESSION ID')
                    .head('RUN CID')
                historyFile.eachLine { table.append(CidHistoryFile.CidRecord.parse(it).toList()) }
                println table.toString()
            } else {
                println("No workflow runs CIDs found.")
            }
        }

        @Override
        void usage() {
            println 'Usage: nextflow cid log'
        }
    }

    class CmdShow implements SubCmd{

        @Override
        String getName() {
            return 'show'
        }

        @Override
        void apply(List<String> args) {
            if (args.size() != 1) {
                println("ERROR: Incorrect number of parameters")
                usage()
                return
            }
            if (!args[0].startsWith(CID_PROT))
                throw new Exception("Identifier is not a CID URL")
            final key = args[0].substring(CID_PROT.size()) + "/$METADATA_FILE"
            final config = new ConfigBuilder()
                    .setOptions(getLauncher().getOptions())
                    .setBaseDir(Paths.get('.'))
                    .build()
            final store = CidStoreFactory.getOrCreate(new Session(config))
            if (store) {
                try {
                    println store.load(key).toString()
                } catch (Throwable e) {
                    println "Error loading ${args[0]}."
                }
            } else {
                println "Error CID store not loaded. Check Nextflow configuration."
            }
        }

        @Override
        void usage() {
            println 'Usage: nextflow cid show <cid reference>'
        }
    }

    class CmdLineage implements SubCmd {

        @Canonical
        class Edge {
            String source
            String destination
            String label
        }

        @Override
        String getName() { 'lineage' }

        @Override
        void apply(List<String> args) {
            if (args.size() != 2) {
                println("ERROR: Incorrect number of parameters")
                usage()
                return
            }
            try {
                final config = new ConfigBuilder()
                    .setOptions(getLauncher().getOptions())
                    .setBaseDir(Paths.get('.'))
                    .build()
                final store = CidStoreFactory.getOrCreate(new Session(config))
                final template = readTemplate()
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
            final key = nodeToRender.substring(CID_PROT.size()) + "/$METADATA_FILE"
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
                    final label = convertToLabel(value.path.toString())
                    lines << "    ${value.path}@{shape: document, label: \"${label}\"}".toString();
                    edges.add(new Edge(value.path.toString(), nodeToRender))
                    return
                }
            }
            final label = convertToLabel(value.toString())
            lines << "    ${value.toString()}@{shape: document, label: \"${label}\"}".toString();
            edges.add(new Edge(value.toString(), nodeToRender))
        }

        protected static String readTemplate() {
            final writer = new StringWriter()
            final res = MermaidHtmlRenderer.class.getResourceAsStream('mermaid.dag.template.html')
            int ch
            while( (ch=res.read()) != -1 ) {
                writer.append(ch as char)
            }
            writer.toString()
        }

        @Override
        void usage() {
            println 'Usage: nextflow cid lineage <cid workflow output > <html output file>'
        }

    }
}
