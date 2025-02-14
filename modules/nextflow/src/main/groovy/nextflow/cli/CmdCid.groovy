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
import nextflow.data.cid.CidStore
import nextflow.data.cid.DefaultCidStore
import nextflow.data.cid.model.DataType
import nextflow.data.config.DataConfig
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins

import java.nio.file.Path
import java.nio.file.Paths

import static nextflow.data.cid.CidObserver.*

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
            final config = new ConfigBuilder()
                    .setOptions(getLauncher().getOptions())
                    .setBaseDir(Paths.get('.'))
                    .build()
            final session = new Session(config)
            final store = session.cidStore
            println store.load("${args[0]}/$METADATA_FILE").toString()
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
                final session = new Session(config)
                final store = session.cidStore
                final template = readTemplate()
                final network = getLineage(store, args[0])
                Path file = Path.of(args[1])
                file.text = template.replace('REPLACE_WITH_NETWORK_DATA', network)
                println("Linage graph for ${args[0]} rendered in ${args[1]}")
            } catch (Throwable e) {
                println("ERROR: rendering lineage graph. ${e.getLocalizedMessage()}")
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
            final slurper = new JsonSlurper()
            final cidObject = slurper.parse(store.load("$nodeToRender/$METADATA_FILE").toString().toCharArray()) as Map
            switch (DataType.valueOf(cidObject.type as String)) {
                case DataType.Output:
                    lines << "    ${nodeToRender}@{shape: document, label: \"${nodeToRender}\"}".toString();
                    final source = cidObject.source as String
                    if (source) {
                        if (source.startsWith(CID_PROT)) {
                            final cid = source.substring(CID_PROT.size())
                            nodes.add(cid)
                            edges.add(new Edge(cid, nodeToRender))
                        } else {
                            lines << "    ${source}@{shape: document, label: \"${source}\"}".toString();
                            edges.add(new Edge(source, nodeToRender))
                        }
                    }

                    break;
                case DataType.WorkflowRun:
                    lines << "${nodeToRender}@{shape: processes, label: \"${cidObject.runName}\"}".toString()
                    final parameters = cidObject.params as Map
                    parameters.values().each {
                        lines << "    ${it}@{shape: document, label: \"${it}\"}".toString();
                        edges.add(new Edge(it.toString(), nodeToRender))
                    }
                    break;
                case DataType.Task:
                    lines << "    ${nodeToRender}@{shape: process, label: \"${cidObject.name}\"}".toString()
                    final parameters = cidObject.inputs as Map<String, String>
                    parameters.values().each { String source ->
                        if (source.startsWith(CID_PROT)) {
                            final cid = source.substring(CID_PROT.size())
                            nodes.add(cid)
                            edges.add(new Edge(cid, nodeToRender))
                        } else {
                            lines << "    ${source}@{shape: document, label: \"${source}\"}".toString();
                            edges.add(new Edge(source, nodeToRender))
                        }
                    }
                    break;
                default:
                    throw new Exception("Unrecognized type reference ${cidObject.type}")
            }
        }

        private String readTemplate() {
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
