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

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cli.CmdBase
import nextflow.pipeline.PipelineManager
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import nextflow.util.StringUtils
import nextflow.util.TestOnly
import nextflow.util.VersionNumber

/**
 * Pipeline list subcommand
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "List installed pipelines")
class CmdPipelineList extends CmdBase {

    @Parameter(names = ['-o', '-output'], description = 'Output mode: text, json')
    String output = 'text'

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
        PipelineCmdSupport.checkOutputMode(output, ['text', 'json'])
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
        final rows = userLevelRows()
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

    /** One row per downloaded commit of each pipeline in the user-level store */
    private static List<Map<String,Object>> userLevelRows() {
        Plugins.init()
        final rows = new ArrayList<Map<String,Object>>()
        for( final String name : AssetManager.list() ) {
            try (final manager = new AssetManager(name)) {
                // a registry name may differ from the Git repository, so report the clone's remote
                final origin = manager.gitConfigRemoteUrl
                final repositoryUrl = origin ? StringUtils.redactUrlPassword(origin) : manager.repositoryUrl
                // an annotated tag and its commit may be downloaded in two clones
                final byCommit = manager.listDownloadedRevisions().groupBy { it.commitId }
                for( final revs : byCommit.values() ) {
                    final rev = revs.find { it.localPath.name == it.commitId } ?: revs[0]
                    final version = revs.find { it.name }?.name
                    rows << row(manager.project, version, rev.commitId, revs.any { !it.clean }, rev.localPath.toPath(), repositoryUrl)
                }
            }
            catch( Exception e ) {
                log.warn "Unable to read pipeline ${name} in the user-level store -- cause: ${e.message ?: e}"
            }
        }
        rows.sort { Map a, Map b -> (a.name as String) <=> (b.name as String) ?: compareVersions(a.version as String, b.version as String) }
        return rows
    }

    /** Compares versions by their numeric parts as numbers and the other parts as text; a revision without a version comes last */
    protected static int compareVersions(String a, String b) {
        if( a == null || b == null )
            return a == b ? 0 : a == null ? 1 : -1
        return new VersionNumber(a) <=> new VersionNumber(b)
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
