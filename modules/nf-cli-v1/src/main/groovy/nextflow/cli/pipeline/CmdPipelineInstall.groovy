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
import groovy.util.logging.Slf4j
import io.seqera.npr.client.RegistryClient
import nextflow.Const
import nextflow.cli.CmdBase
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.module.RegistryClientFactory
import nextflow.pipeline.PipelineGitFetcher
import nextflow.pipeline.PipelineInfo
import nextflow.pipeline.PipelineManager
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import nextflow.script.ScriptFile
import nextflow.util.TestOnly

/**
 * Pipeline install subcommand
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Install a pipeline")
class CmdPipelineInstall extends CmdBase {

    @Parameter(description = "<pipeline>", required = true)
    List<String> args

    @Parameter(names = ["-version"], description = "Pipeline version (default: latest), or the Git branch, tag or commit when installing a repository URL with `-g` (default: the default branch)")
    String version

    // hidden alias of -version, as JCommander can't hide one name of a parameter
    @Parameter(names = ["-r"], hidden = true)
    String revision

    @Parameter(names = ["-force"], description = "Overwrite local modifications (not supported with `-g`)", arity = 0)
    boolean force = false

    @Parameter(names = ['-g', '-global'], description = 'Install the pipeline in the user-level store used by `nextflow run`', arity = 0)
    boolean global

    @TestOnly
    protected RegistryClient client

    @TestOnly
    protected PipelineGitFetcher fetcher

    @TestOnly
    protected Path root

    @Override
    String getName() {
        return 'install'
    }

    @Override
    void run() {
        if( !args || args.size() != 1 )
            throw new AbortOperationException("Unexpected number of parameters")

        version = version ?: revision
        PipelineCmdSupport.checkPipelineName(args[0])
        if( global ) {
            installUserLevel(args[0], version)
            return
        }

        final baseDir = root ?: Path.of('.').toAbsolutePath().normalize()
        final client = this.client ?: RegistryClientFactory.forConfig(PipelineCmdSupport.registryConfig(this, [], baseDir))
        final manager = new PipelineManager(baseDir, client, fetcher ?: new PipelineGitFetcher())
        final result = manager.install(args[0], version, force)
        if( result.alreadyInstalled() )
            println "Pipeline ${result.info().name}@${result.info().version} is already installed"
    }

    private void installUserLevel(String name, String version) {
        if( force )
            throw new AbortOperationException("Option `-force` is not supported with `-global` -- the user-level store keeps each downloaded revision and never overwrites it")
        Plugins.init()
        // the config is read only for a registry address, so that a Git URL never reads it
        RegistryConfig config = null
        final Closure<RegistryConfig> registryConfig = { config ?: (config = PipelineCmdSupport.registryConfig(this)) }
        final isUrl = AssetManager.isUrl(name)
        // a page URL addresses its own registry, a name the configured registries as for project-level installs
        final address = isUrl ? PipelineCmdSupport.registryAddress(name, null, registryConfig) : null
        if( !isUrl || address ) {
            final client = this.client ?: (address ? RegistryClientFactory.forUrls([address.url], address.apiKey) : RegistryClientFactory.forConfig(registryConfig.call()))
            final release = resolveRelease(client, address ? address.name : name, version)
            final label = "${release.name}@${release.version}"
            // a registry release is downloaded at its pinned commit, cached under its registry name as `run` does
            final rev = release.commitId ?: release.revision
            try (final manager = new AssetManager(release.repositoryUrl, rev, Const.DEFAULT_MAIN_FILE_NAME, null, release.name)) {
                final legacyCommit = warnLegacy(manager, label, rev)
                log.info "Checking ${label} ..."
                final upToDate = download(manager, rev, legacyCommit)
                printResult(label, upToDate, manager.getScriptFile())
            }
        }
        else {
            try (final manager = new AssetManager(name)) {
                final legacyCommit = warnLegacy(manager, version ?: 'the latest commit of the checked-out branch', version)
                if( version )
                    manager.setRevision(version)
                log.info "Checking ${manager.getProjectWithRevision()} ..."
                final upToDate = download(manager, version, legacyCommit)
                final script = manager.getScriptFile()
                printResult("${script.projectName}@${version ?: script.revisionInfo.name ?: script.revisionInfo.commitId}", upToDate, script)
            }
        }
    }

    /** Warn that installing into a legacy download replaces its single revision; returns its checked-out commit, or {@code null} when it is not one */
    private static String warnLegacy(AssetManager manager, String label, String revision) {
        if( !manager.isUsingLegacyStrategy() || !manager.localPath?.exists() )
            return null
        final current = manager.getCurrentRevisionAndName()
        if( !revision || (revision != current?.name && !current?.commitId?.startsWith(revision)) )
            log.warn "Pipeline ${manager.project} was downloaded with the legacy layout, which holds a single revision -- installing ${label} replaces the current revision"
        return current?.commitId
    }

    /** Download a revision into the user-level store, returning whether it was already downloaded */
    private static boolean download(AssetManager manager, String revision, String legacyCommit) {
        final result = manager.download(revision)
        manager.updateModules()
        // the result of a legacy download does not tell whether it switched to another commit
        return legacyCommit
            ? manager.getCurrentRevisionAndName()?.commitId == legacyCommit
            : result == 'Already-up-to-date'
    }

    private static void printResult(String label, boolean upToDate, ScriptFile script) {
        if( upToDate )
            println "Pipeline ${label} is already installed"
        else
            println "Pipeline ${label} installed in ${script.localPath}"
    }

    private static PipelineInfo resolveRelease(RegistryClient client, String name, String version) {
        try {
            return PipelineManager.resolveRelease(client, name, version)
        }
        catch( PipelineManager.PipelineNotFoundException e ) {
            throw new AbortOperationException("Pipeline ${name} is not in the registry -- a Git pipeline is installed by its repository URL, e.g. `nextflow pipeline install -g https://github.com/${name}`", e)
        }
    }
}
