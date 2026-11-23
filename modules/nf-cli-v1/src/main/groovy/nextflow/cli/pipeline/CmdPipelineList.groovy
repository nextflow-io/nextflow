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

package nextflow.cli.pipeline

import java.nio.file.Path

import com.beust.jcommander.IParameterValidator
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.beust.jcommander.Parameters
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cli.CmdBase
import nextflow.pipeline.PipelineManager
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import nextflow.util.TestOnly

/**
 * Pipeline list subcommand
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "List installed pipelines")
class CmdPipelineList extends CmdBase {

    @Parameter(names = ['-o', '-output'], description = 'Output mode: text, json', validateWith = OutputModeValidator)
    String output = 'text'

    static class OutputModeValidator implements IParameterValidator {

        private static final List<String> MODES = List.of('text', 'json')

        @Override
        void validate(String name, String value) {
            if( !MODES.contains(value) )
                throw new ParameterException("Output mode must be one of $MODES (found: $value)")
        }
    }

    @Parameter(names = ['-g', '-global'], description = 'List the pipelines in the user-level store used by `nextflow run`', arity = 0)
    boolean global

    @TestOnly
    protected Path root

    @Override
    String getName() {
        return 'list'
    }

    @Override
    void run() {
        if( global ) {
            listUserLevel()
            return
        }
        final manager = new PipelineManager(root ?: Path.of('.').toAbsolutePath().normalize(), null)
        final rows = manager.list().collect { p ->
            row(p.info().name, p.info().version, p.info().commitId, p.info().isModified(p.dir()), p.dir(), p.info().repositoryUrl)
        }
        if( output == 'json' ) {
            printJson(rows)
            return
        }
        if( !rows ) {
            println "No pipelines installed"
            final count = AssetManager.list().size()
            if( count )
                println "${count} ${count == 1 ? 'pipeline is' : 'pipelines are'} installed in the user-level store -- use `nextflow pipeline list -g` to list them"
            return
        }
        printTable(rows)
    }

    private void listUserLevel() {
        Plugins.init()
        final rows = new ArrayList<Map<String,Object>>()
        for( final String name : AssetManager.list() ) {
            try (final manager = new AssetManager(name)) {
                for( final rev : manager.listDownloadedRevisions() )
                    rows << row(manager.project, rev.name, rev.commitId, !rev.clean, rev.localPath.toPath(), manager.repositoryUrl)
            }
            catch( Exception e ) {
                log.warn "Unable to read pipeline ${name} in the user-level store -- cause: ${e.message ?: e}"
            }
        }
        rows.sort { Map a, Map b -> (a.name as String) <=> (b.name as String) ?: compareVersions(a.version as String, b.version as String) }
        if( output == 'json' ) {
            printJson(rows)
            return
        }
        if( !rows ) {
            println "No pipelines installed in the user-level store"
            return
        }
        printTable(rows)
    }

    // a revision without a version comes last
    private static int compareVersions(String a, String b) {
        if( a == b )
            return 0
        if( a == null )
            return 1
        if( b == null )
            return -1
        return a <=> b
    }

    private static Map<String,Object> row(String name, String version, String commit, boolean modified, Path path, String repositoryUrl) {
        return [
            name         : name,
            version      : version,
            commit       : commit,
            status       : modified ? 'modified' : 'ok',
            path         : path?.toString(),
            repositoryUrl: repositoryUrl
        ] as Map<String,Object>
    }

    private static void printJson(List<Map<String,Object>> rows) {
        println JsonOutput.prettyPrint(JsonOutput.toJson(pipelines: rows))
    }

    private static void printTable(List<Map<String,Object>> rows) {
        final table = [['NAME', 'VERSION', 'COMMIT', 'STATUS']]
        rows.each { Map r -> table << [cell(r.name), cell(r.version), cell((r.commit as String)?.take(10)), cell(r.status)] }
        final widths = (0..<3).collect { int i -> table.collect { List<String> it -> it[i].length() }.max() }
        table.each { List<String> cols ->
            println((0..<3).collect { int i -> cols[i].padRight(widths[i] + 2) }.join('') + cols[3])
        }
    }

    // a damaged .pipeline-info may lack any field
    private static String cell(Object value) {
        return value != null ? value.toString() : '-'
    }
}
