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
import io.seqera.npr.client.RegistryException
import nextflow.Const
import nextflow.cli.CmdBase
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
        final config = PipelineCmdSupport.registryConfig(this)
        final isUrl = PipelineCmdSupport.isUrl(name)
        // a page URL addresses its own registry, a name the configured registries as for project-level installs
        final address = isUrl
            ? PipelineCmdSupport.registryAddress(name, null, '-version') { config }
            : new PipelineCmdSupport.RegistryAddress(config.url, name, config.apiKey)
        if( address ) {
            final client = this.client ?: (isUrl ? RegistryClientFactory.forUrls([address.url], address.apiKey) : RegistryClientFactory.forConfig(config))
            final release = resolveRelease(address, version, client)
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
                final legacyCommit = warnLegacy(manager, version ?: 'the default branch', version)
                if( version )
                    manager.setRevision(version)
                log.info "Checking ${manager.getProjectWithRevision()} ..."
                final upToDate = download(manager, version, legacyCommit)
                final script = manager.getScriptFile()
                printResult("${script.projectName}@${version ?: script.revisionInfo.name ?: script.revisionInfo.commitId}", upToDate, script)
            }
        }
    }

    /**
     * Warn that installing a revision of a legacy download replaces its single revision, unless it is the checked-out one
     *
     * @param label The revision as shown to the user
     * @param revision The revision to install: a name or a commit, or {@code null} for the default branch
     * @return the checked-out commit of the legacy download, or {@code null} when it is not one
     */
    private static String warnLegacy(AssetManager manager, String label, String revision) {
        if( !manager.isUsingLegacyStrategy() || !manager.localPath?.exists() )
            return null
        final current = manager.getCurrentRevisionAndName()
        if( !revision || (revision != current?.name && !current?.commitId?.startsWith(revision)) )
            log.warn "Pipeline ${manager.project} was downloaded with the legacy layout, which holds a single revision -- installing ${label} replaces the current revision"
        return current?.commitId
    }

    /**
     * Download a revision into the user-level store
     *
     * @param legacyCommit The checked-out commit of a legacy download, or {@code null}
     * @return whether the revision was already downloaded
     */
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

    private static PipelineInfo resolveRelease(PipelineCmdSupport.RegistryAddress address, String version, RegistryClient client) {
        try {
            return PipelineCmdSupport.resolveRelease(address, version, client)
        }
        catch( AbortOperationException e ) {
            if( isNotFound(e.cause) && (!version || !hasPipeline(client, address.name)) )
                throw new AbortOperationException("Pipeline ${address.name} is not in the registry -- a Git pipeline is installed by its repository URL, e.g. `nextflow pipeline install -g https://github.com/${address.name}`", e)
            throw e
        }
    }

    private static boolean hasPipeline(RegistryClient client, String name) {
        try {
            return client.getPipeline(name) != null
        }
        catch( RegistryException e ) {
            if( isNotFound(e) )
                return false
            throw e
        }
    }

    private static boolean isNotFound(Throwable e) {
        return e instanceof RegistryException && (e as RegistryException).statusCode == 404
    }
}
