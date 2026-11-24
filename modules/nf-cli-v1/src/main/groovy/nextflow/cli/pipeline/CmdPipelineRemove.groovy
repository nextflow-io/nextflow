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
import groovy.transform.CompileStatic
import nextflow.cli.CmdBase
import nextflow.exception.AbortOperationException
import nextflow.pipeline.PipelineManager
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import nextflow.util.TestOnly

/**
 * Pipeline remove subcommand
 */
@CompileStatic
@Parameters(commandDescription = "Remove an installed pipeline")
class CmdPipelineRemove extends CmdBase {

    @Parameter(description = "<pipeline>", required = true)
    List<String> args

    @Parameter(names = ["-version"], description = "Remove only this revision from the user-level store: a version, branch or commit (requires `-g`)")
    String version

    // hidden alias of -version, as JCommander can't hide one name of a parameter
    @Parameter(names = ["-r"], hidden = true)
    String revision

    @Parameter(names = ["-force"], description = "Remove even with local modifications", arity = 0)
    boolean force = false

    @Parameter(names = ['-g', '-global'], description = 'Remove the pipeline from the user-level store used by `nextflow run`', arity = 0)
    boolean global

    @TestOnly
    protected Path root

    @Override
    String getName() {
        return 'remove'
    }

    @Override
    void run() {
        if( !args || args.size() != 1 )
            throw new AbortOperationException("Unexpected number of parameters")
        version = version ?: revision
        PipelineCmdSupport.checkPipelineName(args[0])
        if( global ) {
            removeUserLevel(PipelineCmdSupport.userLevelName(this, args[0]), version)
            return
        }
        if( version )
            throw new AbortOperationException("Option `-version` requires `-global` -- the project-level store holds a single version of each pipeline")
        final manager = new PipelineManager(root ?: Path.of('.').toAbsolutePath().normalize(), null)
        if( manager.remove(args[0], force) ) {
            println "Pipeline ${args[0]} removed"
            return
        }
        println "Pipeline ${args[0]} is not installed"
        if( isInUserLevelStore(args[0]) )
            println "Pipeline ${args[0]} is installed in the user-level store -- use `nextflow pipeline remove -g ${args[0]}` to remove it"
    }

    private static boolean isInUserLevelStore(String name) {
        return AssetManager.list().contains(name)
    }

    private void removeUserLevel(String name, String version) {
        Plugins.init()
        try (final manager = new AssetManager(name)) {
            if( manager.isNotInitialized() ) {
                println "Pipeline ${name} is not installed in the user-level store"
                return
            }
            if( !version ) {
                manager.drop(null, force)
                println "Pipeline ${manager.project} removed"
                return
            }
            final revisions = manager.listDownloadedRevisions()
            final commit = findCommit(manager.project, revisions, version)
            // every clone of the commit, e.g. of an annotated tag and of its commit
            final clones = revisions.findAll { it.commitId == commit }
            if( !force && clones.any { !it.clean } )
                throw new AbortOperationException("Pipeline ${manager.project}@${version} has local modifications -- use -force to remove it")
            if( manager.isUsingLegacyStrategy() ) {
                manager.drop(null, force)
            }
            else {
                for( final rev : clones )
                    manager.drop(rev.localPath.name, force)
                // removing the last clone on disk, readable or not, removes the pipeline
                if( !hasClones(clones[0].localPath.parentFile) )
                    manager.drop(null, force)
            }
            println "Pipeline ${manager.project}@${version} removed"
        }
    }

    private static boolean hasClones(File revisionsDir) {
        return revisionsDir.listFiles()?.any { File it -> it.isDirectory() }
    }

    /** The commit of the downloaded revisions with the given version or branch, else that starts with the given commit or clone name (at least 7 characters) */
    private static String findCommit(String project, List<AssetManager.DownloadedRevision> revisions, String version) {
        final byName = revisions.findAll { it.names.contains(version) }
        final isCommit = version ==~ /[0-9a-f]{7,40}/
        final matches = byName ?: !isCommit ? byName : revisions.findAll { it.commitId.startsWith(version) || it.localPath.name.startsWith(version) }
        if( !matches ) {
            final available = revisions.collect { it.name ?: it.commitId.take(10) }.unique().sort()
            throw new AbortOperationException("Pipeline ${project}@${version} is not installed in the user-level store" + (available ? " -- installed revisions: ${available.join(', ')}" : ''))
        }
        final commits = matches.collect { it.commitId }.unique()
        if( commits.size() > 1 )
            throw new AbortOperationException("Ambiguous revision ${version} for pipeline ${project} -- matches commits: ${commits.collect { it.take(10) }.join(', ')}")
        return commits[0]
    }
}
