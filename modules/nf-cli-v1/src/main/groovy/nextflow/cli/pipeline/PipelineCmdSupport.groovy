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

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.npr.client.RegistryClient
import nextflow.SysEnv
import nextflow.cli.CmdBase
import nextflow.config.ConfigCmdAdapter
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.module.RegistryClientFactory
import nextflow.pipeline.PipelineInfo
import nextflow.pipeline.PipelineManager
import nextflow.scm.AssetManager
import nextflow.scm.HubOptions
import nextflow.scm.ProviderConfig
import nextflow.script.ScriptFile

/**
 * Helpers shared by the pipeline subcommands
 */
@Slf4j
@CompileStatic
class PipelineCmdSupport {

    /**
     * Registry config from the user config, the `nextflow.config` in the base directory, and the given config files
     *
     * @param baseDir The launch directory, the current directory by default
     */
    static RegistryConfig registryConfig(CmdBase cmd, List<Path> configFiles = [], Path baseDir = null) {
        final config = new ConfigCmdAdapter()
            .setOptions(cmd.launcher.options)
            .setUserConfigFiles(configFiles)
            .setBaseDir(baseDir ?: Path.of('.').toAbsolutePath().normalize())
            .build()
        return config.navigate('registry') as RegistryConfig ?: new RegistryConfig()
    }

    /**
     * The name of a pipeline in the user-level store, where a registry pipeline is cached under its registry name.
     * Used by commands that take no version.
     *
     * @param name A pipeline name or a registry page URL
     */
    static String userLevelName(CmdBase cmd, String name) {
        checkPipelineName(name, null)
        final address = registryAddress(name, null, null) { registryConfig(cmd) }
        return address ? address.name : name
    }

    /**
     * Fails for an SSH Git URL, which is not supported, and for a pipeline name with an {@code @version} suffix
     *
     * @param versionOption The option that gives the version, suggested instead of the suffix, or {@code null} when the command takes no version
     */
    static void checkPipelineName(String name, String versionOption = '-version') {
        final ssh = name =~ /^(?:ssh:\/\/)?(?:[\w.-]+@)?([\w.-]+)(?::\d+\/|:|\/)(.+?)(?:\.git)?\/?$/
        if( name.startsWith('ssh://') || (name =~ /^[\w.-]+@[\w.-]+:/).find() ) {
            final hint = ssh.matches() ? " -- use the https URL instead, e.g. `https://${ssh.group(1)}/${ssh.group(2)}`" : ' -- use the https URL instead'
            throw new AbortOperationException("SSH Git URLs are not supported: ${name}${hint}")
        }
        if( !isUrl(name) && name.contains('@') ) {
            final at = name.lastIndexOf('@')
            final hint = versionOption
                ? "give the version with the `${versionOption}` option, e.g. `${name.substring(0, at)} ${versionOption} ${name.substring(at + 1)}`"
                : "give the pipeline as `scope/name`, e.g. `${name.substring(0, at)}`"
            throw new AbortOperationException("Not a valid pipeline name: ${name} -- ${hint}")
        }
    }

    /** Whether a pipeline is given by a URL rather than by name, as {@link AssetManager} tells them apart */
    static boolean isUrl(String name) {
        return name.startsWith('http://') || name.startsWith('https://') || name.startsWith('file:/')
    }

    /** A registry pipeline address: the registry API URL, the pipeline name, and the API key to send */
    @Canonical
    static class RegistryAddress {
        String url
        String name
        String apiKey
    }

    /**
     * The registry addressed by a registry page URL, or by a registry host as hub ({@code -hub} or {@code NXF_HUB}).
     * A hub that is a Git provider name keeps addressing Git, as before.
     *
     * @return the registry address, or {@code null} when {@code name} and {@code hub} address Git
     */
    static RegistryAddress registryAddress(String name, String hub, Closure<RegistryConfig> registryConfig) {
        return registryAddress(name, hub, '-r', registryConfig)
    }

    /**
     * @param versionOption The option that gives the version, suggested when the address includes one, or {@code null} when the command takes no version
     */
    static RegistryAddress registryAddress(String name, String hub, String versionOption, Closure<RegistryConfig> registryConfig) {
        final isUrl = name.startsWith('http://') || name.startsWith('https://')
        hub = hub ?: SysEnv.get('NXF_HUB')
        // cheap checks first, so that Git addresses never read the config
        if( !isUrl && !hub )
            return null
        final uri = isUrl ? toUri(name) : null
        // a registry URL is a pipeline page
        if( isUrl && !uri?.path?.contains('/pipelines/') )
            return null
        // a Git provider host is known without the launch config, e.g. a GitLab group named `pipelines`
        if( isUrl && isGitProviderHost(uri) )
            return null
        final providers = isUrl ? null : ProviderConfig.getDefaultProviderNames()
        if( !isUrl && hub in providers )
            return null
        final config = registryConfig.call()
        // the default public registry is always recognized, even when not configured
        final registries = config.allUrls + [RegistryConfig.DEFAULT_REGISTRY_URL]
        // a URL and a hub name a registry by its host, or host and port
        final api = registryOfHost(isUrl ? uri.authority : hub, registries)
        // an address names one registry, and the API key belongs to the primary one
        final apiKey = api && api == config.url ? config.apiKey : null
        if( isUrl ) {
            if( !api )
                return null
            // a registry URL must be a pipeline page; the scheme is ignored since requests use the configured API URL
            final prefix = (URI.create(pageUrl(api)).path ?: '') + '/pipelines/'
            if( !uri.path.startsWith(prefix) )
                return null
            return newAddress(api, uri.path.substring(prefix.length()).replaceFirst(/\/$/, ''), apiKey, versionOption)
        }
        if( !api ) {
            final hosts = registries.collect { URI.create(it).authority }.unique()
            throw new AbortOperationException("Unknown repository provider or registry: `$hub` -- providers: ${providers.join(', ')}; registries: ${hosts.join(', ')}")
        }
        return newAddress(api, name, apiKey, versionOption)
    }

    // checked here rather than when resolving the release, so that it applies offline too
    private static RegistryAddress newAddress(String api, String name, String apiKey, String versionOption) {
        if( name.contains('@') )
            throw new AbortOperationException("Not a valid registry pipeline address: ${name}" + (versionOption ? " -- give the version with the `${versionOption}` option" : ''))
        return new RegistryAddress(api, name, apiKey)
    }

    private static boolean isGitProviderHost(URI uri) {
        final domains = ProviderConfig.getDefaultProviderDomains()
        return domains.any { String it -> it.equalsIgnoreCase(uri.host) || it.equalsIgnoreCase(uri.authority) }
    }

    /**
     * Resolve the release of a registry address, only against the registry it names.
     *
     * @param revision The registry version, or {@code null} for the latest release
     */
    static PipelineInfo resolveRelease(RegistryAddress address, String revision, RegistryClient client = null) {
        final release = PipelineManager.resolveRelease(client ?: RegistryClientFactory.forUrls([address.url], address.apiKey), address.name, revision)
        log.info "Using registry release ${address.name}@${release.version}" + (release.commitId ? " (commit ${release.commitId.take(10)})" : '')
        return release
    }

    /** A revision downloaded into the user-level store */
    @Canonical
    static class UserLevelDownload {
        /** The download result, e.g. {@code Already-up-to-date} when the revision was already downloaded */
        String message
        ScriptFile script

        boolean isUpToDate() {
            return message == 'Already-up-to-date'
        }
    }

    /**
     * Download a registry release into the user-level store at its pinned commit, cached under its registry name as `run` does
     */
    static UserLevelDownload downloadRelease(PipelineInfo release, HubOptions hubOpts, Integer deep = null, boolean install = false) {
        final rev = release.commitId ?: release.revision
        try (final manager = new AssetManager(release.repositoryUrl, rev, null, hubOpts, release.name)) {
            final legacyCommit = install ? warnLegacy(manager, "${release.name}@${release.version}", rev) : null
            log.info "Checking ${release.name}@${release.version} ..."
            final result = legacySwitch(manager, legacyCommit, manager.download(rev, deep))
            manager.updateModules()
            return new UserLevelDownload(result, manager.getScriptFile(release.workflowPath))
        }
    }

    /**
     * Download a Git revision of a project into the user-level store
     *
     * @param project A project name or a Git repository URL
     * @param revision A branch, tag or commit, or {@code null} for the default branch
     * @param migrate Whether to migrate a legacy download to the multi-revision strategy
     * @param install Whether the download is an install, whose warning for a legacy download does not suggest {@code -migrate}
     */
    static UserLevelDownload downloadRevision(String project, String revision, HubOptions hubOpts, Integer deep = null, boolean migrate = false, boolean install = false) {
        try (final manager = new AssetManager(project, hubOpts)) {
            final legacyCommit = install ? warnLegacy(manager, revision ?: 'the default branch', revision) : null
            if( !install && manager.isUsingLegacyStrategy() ) {
                if( migrate ) {
                    log.info "Migrating ${project} revision ${revision} to multi-revision strategy"
                    manager.setStrategyType(AssetManager.RepositoryStrategyType.MULTI_REVISION)
                } else {
                    log.warn "The local asset for ${project} does not support multi-revision - Pulling with legacy strategy\n" +
                        "Consider updating the project ${project} using '-migrate' option"
                }
            }

            if( revision )
                manager.setRevision(revision)

            log.info "Checking ${manager.getProjectWithRevision()} ..."

            final result = legacySwitch(manager, legacyCommit, manager.download(revision, deep))
            manager.updateModules()
            return new UserLevelDownload(result, manager.getScriptFile())
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

    // the result of a legacy download does not tell whether it switched to another commit
    private static String legacySwitch(AssetManager manager, String before, String result) {
        if( !before )
            return result
        final after = manager.getCurrentRevisionAndName()?.commitId
        return after == before ? 'Already-up-to-date' : "checked out at ${after}"
    }

    /** Base URL of the registry web pages: the API URL without {@code /api} */
    static String pageUrl(String apiUrl) {
        return apiUrl.replaceFirst(/\/api\/?$/, '')
    }

    private static String registryOfHost(String host, Collection<String> registryUrls) {
        return registryUrls.find { final uri = URI.create(it); uri.host?.equalsIgnoreCase(host) || uri.authority?.equalsIgnoreCase(host) }
    }

    /** A URL that {@link URI} can't parse is not a registry URL, so it is left to Git */
    private static URI toUri(String url) {
        try {
            return URI.create(url)
        }
        catch( IllegalArgumentException e ) {
            return null
        }
    }
}
