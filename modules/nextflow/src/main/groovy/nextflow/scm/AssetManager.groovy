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

package nextflow.scm

import static nextflow.Const.*
import static GitReferenceHelper.*

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import nextflow.cli.HubOptions
import nextflow.config.Manifest
import nextflow.config.ConfigParserFactory
import nextflow.exception.AbortOperationException
import nextflow.exception.AmbiguousPipelineNameException
import nextflow.script.ScriptFile
import nextflow.SysEnv
import nextflow.util.IniFile
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.api.errors.RefNotFoundException
import org.eclipse.jgit.lib.Ref
import org.eclipse.jgit.lib.Repository
import org.eclipse.jgit.merge.MergeStrategy
/**
 * Handles operation on remote and local installed pipelines.
 * Uses the strategy pattern to support different ways to manage local installations.
 * Current available {@link RepositoryStrategy}:
 * - {@link LegacyRepositoryStrategy}: This is the traditional approach where each project gets a full git clone.
 * - {@Link MultiRevisionRepositoryStrategy}: This approach allows multiple revisions to coexist efficiently by sharing objects
 * through a bare repository and creating lightweight clones for each commit.
 *
 * A {@link AssetManager.RepositoryStatus} is defined according to the status of the project folder (`localRootPath`).
 * It is used to automatically select the {@link RepositoryStrategy}. The {@link LegacyRepositoryStrategy} will be selected for LEGACY_ONLY status,
 * and the @Link MultiRevisionRepositoryStrategy} for other statuses (UNINNITIALIZED, BARE_ONLY and HYBRID)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
class AssetManager implements Closeable {

    /**
     * The folder all pipelines scripts are installed
     */
    @PackageScope
    static File root = DEFAULT_ROOT

    /**
     * The pipeline name. It must be in the form {@code username/repo} where 'username'
     * is a valid user name or organization account, while 'repo' is the repository name
     * containing the pipeline code
     */
    private String project

    private Git _git

    private String mainScript

    private RepositoryProvider provider

    private String hub

    private List<ProviderConfig> providerConfigs

    /**
     * The repository strategy (legacy or multi-revision)
     */
    private RepositoryStrategy strategy

    /**
     * Create a new asset manager object with default parameters
     */
    AssetManager() {
        this.providerConfigs = ProviderConfig.createDefault()
    }

    /**
     * Create a new asset manager with the specified pipeline name
     *
     * @param pipeline The pipeline to be managed by this manager e.g. {@code nextflow-io/hello}
     */
    AssetManager( String pipelineName, HubOptions cliOpts = null) {
        assert pipelineName
        // read the default config file (if available)
        def config = ProviderConfig.getDefault()
        // build the object
        build(pipelineName, config, cliOpts)
    }

    AssetManager( String pipelineName, Map config ) {
        assert pipelineName
        // build the object
        build(pipelineName, config)
    }

    AssetManager( String pipelineName, String revision, HubOptions cliOpts = null ) {
        assert pipelineName
        // build the object
        def config = ProviderConfig.getDefault()
        // build the object
        build(pipelineName, config, cliOpts, revision)
    }

    /**
     * Build the asset manager internal data structure
     *
     * @param pipelineName A project name or a project repository Git URL
     * @param config A {@link Map} holding the configuration properties defined in the {@link ProviderConfig#DEFAULT_SCM_FILE} file
     * @param cliOpts User credentials provided on the command line. See {@link HubOptions} trait
     * @return The {@link AssetManager} object itself
     */
    @PackageScope
    AssetManager build( String pipelineName, Map config = null, HubOptions cliOpts = null, String revision = null ) {

        this.providerConfigs = ProviderConfig.createFromMap(config)

        this.project = resolveName(pipelineName)

        if( !isValidProjectName(this.project) ) {
            throw new IllegalArgumentException("Not a valid project name: ${this.project}")
        }
        // Initialize strategy based on environment and repository state
        initStrategy(revision)
        this.hub = checkHubProvider(cliOpts)
        this.provider = createHubProvider(hub)

        if( revision ){
            setRevision(revision)
        }

        strategy.setProvider(this.provider)

        setupCredentials(cliOpts)

        validateProjectDir()

        return this
    }

    /**
     * Update the repository strategy based on environment and current state, or create a new one if not exists
     */
    private void updateStrategy(RepositoryStrategyType type = null) {
        def revision = getRevision()
        if( strategy )
            strategy.close()
        if( !type )
            type = selectStrategyType()
        if( type == RepositoryStrategyType.LEGACY )
            strategy = new LegacyRepositoryStrategy(project)
        else if( type == RepositoryStrategyType.MULTI_REVISION )
            strategy = new MultiRevisionRepositoryStrategy(project, revision)
        else
            throw new AbortOperationException("Not supported strategy $type")
        if( revision )
            setRevision(revision)
        strategy.setProvider(provider)
    }

    private void initStrategy(String revision = null) {
        def type = selectStrategyType()
        if( type == RepositoryStrategyType.LEGACY )
            strategy = new LegacyRepositoryStrategy(project)
        else if( type == RepositoryStrategyType.MULTI_REVISION )
            strategy = new MultiRevisionRepositoryStrategy(project, revision)
        else
            throw new AbortOperationException("Not supported strategy $type")
        if( revision )
            setRevision(revision)
        strategy.setProvider(provider)
    }

    private RepositoryStrategyType selectStrategyType() {
        // Check environment variable for legacy mode
        if( SysEnv.get('NXF_SCM_LEGACY') as boolean ) {
            log.warn "Forcing to use legacy repository strategy (NXF_SCM_LEGACY is set to true)"
            return RepositoryStrategyType.LEGACY
        }

        if (isOnlyLegacy()){
            log.debug "Legacy repository detected, selecting legacy strategy"
            return RepositoryStrategyType.LEGACY
        } else {
            log.debug "Selecting multi-revision strategy"
            return RepositoryStrategyType.MULTI_REVISION
        }
    }

    /**
     * Set the repository strategy type explicitly
     *
     * @param useLegacy If true, use legacy strategy; otherwise use multi-revision
     */
    void setStrategyType(RepositoryStrategyType type) {
        if( (type == RepositoryStrategyType.LEGACY && isUsingLegacyStrategy()) ||
            (type == RepositoryStrategyType.MULTI_REVISION && isUsingMultiRevisionStrategy()) ) {
            //nothing to do
            return
        }
        if( type == RepositoryStrategyType.LEGACY && isUsingMultiRevisionStrategy() ) {
            log.warn1 "Switching to legacy SCM repository management"
            updateStrategy(RepositoryStrategyType.LEGACY)
        } else if( type == RepositoryStrategyType.MULTI_REVISION && SysEnv.get('NXF_SCM_LEGACY') ) {
            log.warn1 "Received a request to change to multi-revision SCM repository management and NXF_SCM_LEGACY is defined. Keeping legacy strategy"
            updateStrategy(RepositoryStrategyType.LEGACY)
        } else {
            log.debug "Switching to multi-revision repository strategy"
            updateStrategy(RepositoryStrategyType.MULTI_REVISION)
        }
    }

    /**
     * Check if using legacy strategy
     */
    boolean isUsingLegacyStrategy() {
        return strategy instanceof LegacyRepositoryStrategy
    }

    /**
     * Check if using multi-revision strategy
     */
    boolean isUsingMultiRevisionStrategy() {
        return strategy instanceof MultiRevisionRepositoryStrategy
    }

    /**
     * Check if a multi-revision or legacy local asset exist.
     * @return 'true' if none of them exist, otherwise 'false'.
     */
    boolean isNotInitialized() {
        final hasMultiRevision = MultiRevisionRepositoryStrategy.checkProject(root, project)
        final hasLegacy = LegacyRepositoryStrategy.checkProject(root, project)

        return !hasMultiRevision && !hasLegacy
    }

    /**
     * Check if repository has only a local legacy asset
     * @return 'true' if legacy asset exists and multi-revision doesn't exist, otherwise 'false'.
     */
    boolean isOnlyLegacy(){
        final hasMultiRevision = MultiRevisionRepositoryStrategy.checkProject(root, project)
        final hasLegacy = LegacyRepositoryStrategy.checkProject(root, project)

        return hasLegacy && !hasMultiRevision
    }


    /**
     * Set the revision for multi-revision strategy
     */
    AssetManager setRevision(String revision) {
        if( provider )
            provider.setRevision(revision)
        if( revision && strategy instanceof MultiRevisionRepositoryStrategy ) {
            ((MultiRevisionRepositoryStrategy) strategy).setRevision(revision)
        }
        return this
    }

    /**
     * Get the current revision (multi-revision strategy only)
     */
    private String getRevision() {
        if( provider )
            return provider.getRevision()
        if( strategy && strategy instanceof MultiRevisionRepositoryStrategy ) {
            return ((MultiRevisionRepositoryStrategy) strategy).getRevision()
        }
        return null
    }

    /**
     * Get the project name with revision appended (for multi-revision strategy)
     */
    String getProjectWithRevision() {
        def rev = getRevision()
        return project + (rev ? ':' + rev : '')
    }

    @PackageScope
    File getLocalGitConfig() {
        return strategy?.getLocalGitConfig() ?: null
    }

    @PackageScope
    AssetManager setProject(String name) {
        this.project = name
        if( !strategy )
            initStrategy(getRevision())
        strategy.setProject(name)
        return this
    }

    protected RepositoryProvider getProvider() {
        return provider
    }

    /**
     * Sets the user credentials on the {@link RepositoryProvider} object
     *
     * @param cliOpts The user credentials specified on the program command line. See {@code HubOptions}
     */
    @PackageScope
    void setupCredentials( HubOptions cliOpts ) {
        if( cliOpts?.hubUser ) {
            cliOpts.hubProvider = hub
            final user = cliOpts.getHubUser()
            final pwd = cliOpts.getHubPassword()
            provider.setCredentials(user, pwd)
        }
    }


    @PackageScope
    boolean isValidProjectName( String projectName ) {
        projectName =~~ /.+\/.+/
    }

    /**
     * Verifies that the project hub provider eventually specified by the user using the {@code -hub} command
     * line option or implicitly by entering a repository URL, matches with clone URL of a project already cloned (downloaded).
     */
    @PackageScope
    void validateProjectDir() {

        if( !localPath.exists() ) {
            return
        }

        // if project dir exists it must contain the Git config file
        final configProvider = guessHubProviderFromGitConfig(true)
        if( !configProvider )
            throw new IllegalStateException("Cannot find a provider config for repository at path: $localPath")

        // checks that the selected hub matches with the one defined in the git config file
        if( hub != configProvider ) {
            throw new AbortOperationException("A project with name: `$localPath` has already been downloaded from a different provider: `$configProvider`")
        }

    }

    /**
     * Find out the "hub provider" (i.e. the platform on which the remote repository is stored
     * for example: github, bitbucket, etc) and verifies that it is a known provider.
     *
     * @param cliOpts The user hub info provider as command line options. See {@link HubOptions}
     * @return The name of hub name e.g. {@code github}, {@code bitbucket}, etc.
     */
    @PackageScope
    String checkHubProvider( HubOptions cliOpts ) {

        def result = hub
        if( !result )
            result = cliOpts?.getHubProvider()
        if( !result )
            result = guessHubProviderFromGitConfig()
        if( !result )
            result = DEFAULT_HUB

        def providerNames = providerConfigs.collect { it.name }
        if( !providerNames.contains(result)) {
            def matches = providerNames.closest(result) ?: providerNames
            def message = "Unknown repository provider: `$result`'. Did you mean?\n" + matches.collect { "  $it"}.join('\n')
            throw new AbortOperationException(message)
        }

        return result
    }

    /**
     * Given a project name or a repository URL returns a fully qualified project name.
     *
     * @param name A project name or URL e.g. {@code cbcrg/foo} or {@code https://github.com/cbcrg/foo.git}
     * @return The fully qualified project name e.g. {@code cbcrg/foo}
     */
    @PackageScope
    String resolveName( String name ) {
        assert name

        //
        // check if it's a repository fully qualified URL e.g. https://github.com/foo/bar
        //
        def project = resolveNameFromGitUrl(name)
        if( project )
            return project

        //
        // otherwise it must be a canonical repository name e.g. user/project
        //
        if( ['./','../', '/' ].any(it->name.startsWith(it)) )
            throw new AbortOperationException("Not a valid project name: $name")

        def parts = name.split('/') as List<String>
        def last = parts[-1]
        if( last.endsWith('.nf') || last.endsWith('.nxf') ) {
            if( parts.size()==1 )
                throw new AbortOperationException("Not a valid project name: $name")

            if( parts.size()==2 ) {
                mainScript = last
                parts = [ parts.first() ]
            }
            else {
                mainScript = parts[2..-1].join('/')
                parts = parts[0..1]
            }
        }

        if( parts.size() == 2 ) {
            return parts.join('/')
        }
        else if( parts.size()>2 ) {
            throw new AbortOperationException("Not a valid project name: $name")
        }
        else {
            name = parts[0]
        }

        def qualifiedName = find(name)
        if( !qualifiedName ) {
            return "$DEFAULT_ORGANIZATION/$name".toString()
        }

        if( qualifiedName instanceof List ) {
            final msg = "Which one do you mean?\n${qualifiedName.join('\n')}"
            throw new AmbiguousPipelineNameException(msg, qualifiedName)
        }

        return qualifiedName
    }

    String getProject() { project }

    String getHub() { hub }

    @PackageScope
    String resolveNameFromGitUrl( String repository ) {

        final isUrl = repository.startsWith('http://') || repository.startsWith('https://') || repository.startsWith('file:/')
        if( !isUrl )
            return null

        try {
            def url = new GitUrl(repository)

            def result
            if( url.protocol == 'file' ) {
                this.hub = "file:${url.domain}"
                providerConfigs << new ProviderConfig(this.hub, [path:url.domain])
                result = "local/${url.path}"
            }
            else {
                // find the provider config for this server
                final config = RepositoryFactory.getProviderConfig(providerConfigs, url)
                if( config ) {
                    if( !providerConfigs.contains(config) )
                        providerConfigs.add(config)
                    this.hub = config.name
                    result = config.resolveProjectName(url.path)
                }
                else {
                    result = url.path.stripStart('/')
                }
            }
            log.debug "Repository URL: $repository; Project: $result; Hub provider: $hub"
            return result
        }
        catch( IllegalArgumentException e ) {
            log.debug "Cannot parse Git URL: $repository -- cause: ${e.message}"
        }
        return null
    }


    /**
     * Creates the RepositoryProvider instance i.e. the object that manages the interaction with
     * the remote SCM server (e.g. GitHub, GitLab, etc) using the platform provided API
     *
     * @param providerName The name of the provider e.g. {@code github}, {@code gitlab}, etc
     * @return
     */
    @PackageScope
    RepositoryProvider createHubProvider(String providerName) {

        final config = providerConfigs.find { it.name == providerName }
        if( !config )
            throw new AbortOperationException("Unknown repository configuration provider: $providerName")

        return RepositoryFactory.newRepositoryProvider(config, project)
    }

    AssetManager setLocalPath(File path) {
        if (!strategy)
            initStrategy(getRevision())
        strategy.setLocalPath(path)
        return this
    }

    AssetManager checkValidRemoteRepo(String revision=null) {
        // Configure the git provider to use the required revision as source for all needed remote resources:
        // - config if present in repo (nextflow.config by default)
        // - main script (main.nf by default)
        provider.revision = revision
        final scriptName = getMainScriptName()
        provider.validateFor(scriptName)
        return this
    }

    String getGitRepositoryUrl() {
        return strategy?.getGitRepositoryUrl() ?: provider.getCloneUrl()
    }

    File getLocalPath() {
        return strategy?.getLocalPath()
    }

    File getProjectPath() {
        return strategy?.getProjectPath()
    }

    ScriptFile getScriptFile(String scriptName=null) {

        def result = new ScriptFile(getMainScriptFile(scriptName))
        result.revisionInfo = getCurrentRevisionAndName()
        result.repository = getRepositoryUrl()
        result.localPath = localPath.toPath()
        result.projectName = project

        return result
    }

    File getMainScriptFile(String scriptName=null) {
        if( !localPath.exists() ) {
            throw new AbortOperationException("Unknown project folder: $localPath")
        }

        def mainScript = scriptName ?: getMainScriptName()
        def result = new File(localPath, mainScript)
        if( !result.exists() )
            throw new AbortOperationException("Missing project main script: $result")

        return result
    }

    String getMainScriptName() {
        return mainScript ?: getManifest().getMainScript()
    }

    String getHomePage() {
        getManifest().getHomePage() ?: provider.getRepositoryUrl()
    }

    String getRepositoryUrl() {
        provider?.getRepositoryUrl()
    }

    String getDefaultBranch() {
        // if specified in manifest, that takes priority
        // otherwise look for a symbolic ref (refs/remotes/origin/HEAD)
        return getManifest().getDefaultBranch()
                ?: strategy?.getRemoteDefaultBranch()
                ?: DEFAULT_BRANCH
    }

    @Memoized
    Manifest getManifest() {
        getManifest0()
    }


    protected Manifest getManifest0() {
        String text = null
        ConfigObject result = null
        try {
            text = localPath.exists() ? new File(localPath, MANIFEST_FILE_NAME).text : provider.readText(MANIFEST_FILE_NAME)
        }
        catch( FileNotFoundException e ) {
            log.debug "Project manifest does not exist: ${e.message}"
        }
        catch( Exception e ) {
            log.warn "Cannot read project manifest -- Cause: ${e.message ?: e}", e
        }

        if( text ) try {
            def config = ConfigParserFactory.create().setIgnoreIncludes(true).setStrict(false).parse(text)
            result = (ConfigObject)config.manifest
        }
        catch( Exception e ) {
            log.warn "Cannot read project manifest -- Cause:  ${e.message ?: e}"
        }

        // by default return an empty object
        if( result == null )
            result = new ConfigObject()

        return new Manifest(result)
    }

    Path getConfigFile() {
        if( localPath.exists() ) {
            return new File(localPath, MANIFEST_FILE_NAME).toPath()
        }
        else {
            try {
                // try to read the config file
                provider.readBytes(MANIFEST_FILE_NAME)
                // no error => exist, return a path for it
                return new ProviderPath(provider, MANIFEST_FILE_NAME)
            }
            catch (Exception e) {
                provider.validateRepo()
                log.debug "Cannot retrieve remote config file -- likely it does not exist"
                return null
            }
        }
    }


    String getBaseName() {
        def result = project.tokenize('/')
        if( result.size() > 2 ) throw new IllegalArgumentException("Not a valid project name: $project")
        return result.size()==1 ? result[0] : result[1]
    }

    boolean isLocal() {
        return localPath && localPath.exists()
    }

    /**
     * @return {@code true} when the local project path exists and contains at least the default script
     *      file (i.e. main.nf) or the nextflow manifest file (i.e. nextflow.config)
     */
    boolean isRunnable() {
        isLocal() && ( new File(localPath,DEFAULT_MAIN_FILE_NAME).exists() || new File(localPath,MANIFEST_FILE_NAME).exists() )
    }

    /**
     * @return true if no differences exist between the working-tree, the index,
     *         and the current HEAD, false if differences do exist
     */
    boolean isClean() {
        return strategy?.isClean() ?: true
    }

    /**
     * Close the underlying Git repository
     */
    void close() {
        if( _git ) {
            _git.close()
            _git = null
        }
        if( strategy ) {
            strategy.close()
        }
    }

    /**
     * @return The list of available pipelines
     */
    static List<String> list() {
        log.debug "Listing projects in folder: $root"

        if( !root.exists() )
            return []

        def result = new HashSet<String>()
        appendProjectsFromDir(root, result)
        return result.toList().sort()
    }

    /**
     * Append found projects to results considering the new multi-revision projects in '.repos' subdirectory
     * @param dir
     * @param result
     */
    private static void appendProjectsFromDir(File dir, HashSet result) {
        dir.eachDir { File org ->
            if( org.getName() == MultiRevisionRepositoryStrategy.REPOS_SUBDIR ) {
                appendProjectsFromDir(org, result)
            } else {
                org.eachDir { File it ->
                    log.debug("Adding ${org.getName()}/${it.getName()}")
                    result << "${org.getName()}/${it.getName()}".toString()
                }
            }
        }
    }

    static protected def find( String name ) {
        def exact = []
        def partial = []

        list().each {
            def items = it.split('/')
            if( items[1] == name )
                exact << it
            else if( items[1].startsWith(name ) )
                partial << it
        }

        def list = exact ?: partial
        return list.size() ==1 ? list[0] : list
    }


    protected Git getGit() {
        if( strategy ) {
            return strategy.getGit()
        }
        // Fallback to legacy behavior if strategy not initialized
        if( !_git ) {
            _git = Git.open(localPath)
        }
        return _git
    }

    /**
     * Download a pipeline from a remote Github repository
     *
     * @param revision The revision to download
     * @param deep Optional depth for shallow clones
     * @return A message representing the operation result
     */
    String download(String revision = null, Integer deep = null) {
        assert project
        if( !strategy ) {
            throw new IllegalStateException("Strategy not initialized")
        }
        // If it is a new download check is a valid repository
        if( !localPath.exists() ) {
            checkValidRemoteRepo(revision)
        }
        return strategy.download(revision, deep, manifest)
    }

    /**
     * Clone a pipeline from a remote pipeline repository to the specified folder
     *
     * @param directory The folder when the pipeline will be cloned
     * @param revision The revision to be cloned. It can be a branch, tag, or git revision number
     */
    void clone(File directory, String revision = null, Integer deep=null) {

        def clone = Git.cloneRepository()
        def uri = getGitRepositoryUrl()
        log.debug "Clone project `$project` -- Using remote URI: ${uri} into: $directory"

        if( !uri )
            throw new AbortOperationException("Cannot find the specified project: $project")

        clone.setURI(uri)
        clone.setDirectory(directory)
        clone.setDepth(1)
        clone.setCloneSubmodules(manifest.recurseSubmodules)
        if( provider.hasCredentials() )
            clone.setCredentialsProvider( provider.getGitCredentials() )

        if( revision )
            clone.setBranch(revision)
        if( deep )
            clone.setDepth(deep)
        clone.call().close()
    }

    /**
     * @return The symbolic name of the current revision i.e. the current checked out branch or tag
     */
    String getCurrentRevision() {
        if( !strategy ) {
            throw new IllegalStateException("Strategy not initialized")
        }
        return strategy.getCurrentRevision()

    }

    RevisionInfo getCurrentRevisionAndName() {
        if( !strategy ) {
            throw new IllegalStateException("Strategy not initialized")
        }
        return strategy.getCurrentRevisionAndName()
    }

    /**
     * @return A list of existing branches and tags names. For example
     * <pre>
     *     * master (default)
     *       patch-x
     *       v1.0 (t)
     *       v1.1 (t)
     * </pre>
     *
     * The star character on the left highlight the current revision, the string {@code (default)}
     *  ticks that it is the default working branch (usually the master branch), while the string {@code (t)}
     *  shows that the revision is a git tag (instead of a branch)
     */
    @Deprecated
    List<String> getRevisions(int level) {

        def pulled = strategy.listDownloadedCommits()
        def master = getDefaultBranch()

        List<String> branches = getBranchList()
            .findAll { it.name.startsWith('refs/heads/') || isRemoteBranch(it) }
            .unique { shortenRefName(it.name) }
            .collect { Ref it -> refToString(it, pulled, master, false, level) }

        List<String> tags = getTagList()
                .findAll  { it.name.startsWith('refs/tags/') }
                .collect { refToString(strategy.peel(it) , pulled, master, true, level) }

        def result = new ArrayList<String>(branches.size() + tags.size())
        result.addAll(branches)
        result.addAll(tags)
        return result
    }

    List<RevisionInfo> getRemoteRevisions() {
        final result = new ArrayList(50)
        for( def branch : provider.getBranches() ) {
            result.add(new RevisionInfo(branch.commitId, branch.name, RevisionInfo.Type.BRANCH))
        }
        for( def tag : provider.getTags() ) {
            result.add(new RevisionInfo(tag.commitId, tag.name, RevisionInfo.Type.TAG))
        }
        return result
    }

    Map getBranchesAndTags(boolean checkForUpdates) {
        final result = [:]
        final currentCommits = strategy.listDownloadedCommits()
        final master = getDefaultBranch()

        final branches = []
        final tags = []
        final pulled = new LinkedList<Map>()
        log.debug("Current commits : $currentCommits")
        Map<String, Ref> remote = checkForUpdates ? strategy.lsRemote(false) : null
        getBranchList()
            .findAll { it.name.startsWith('refs/heads/') || isRemoteBranch(it) }
            .unique { shortenRefName(it.name) }
            .each {
                final map = refToMap(it, remote)
                if( isRefInCommits(it, currentCommits) )
                    pulled << map
                branches << map
            }

        remote = checkForUpdates ? strategy.lsRemote(true) : null
        getTagList()
            .findAll { it.name.startsWith('refs/tags/') }
            .each {
                Ref ref = strategy.peel(it)
                final map = refToMap(ref, remote)
                if( isRefInCommits(ref, currentCommits) )
                    pulled << map
                tags << map
            }

        result.pulled = pulled.collect { it.name }  // current pulled revisions
        result.master = master      // master branch name
        result.branches = branches  // collection of branches
        result.tags = tags          // collect of tags
        return result
    }

    @Memoized
    protected List<Ref> getBranchList() {
        return strategy?.getBranchList() ?: []
    }

    @Memoized
    protected List<Ref> getTagList() {
        return strategy?.getTagList() ?: []
    }

    protected String refToString(Ref ref, List<String> downloaded, String master, boolean tag, int level ) {

        def result = new StringBuilder()
        def name = shortenRefName(ref.name)
        def coChar = isOnlyLegacy() ? '*' : '>'
        result << ( isRefInCommits(ref, downloaded) ? coChar : ' ')

        if( level ) {
            def obj = ref.getPeeledObjectId() ?: ref.getObjectId()
            result << ' '
            result << formatObjectId(obj, level == 1)
        }

        result << ' ' << name

        if( tag )
            result << ' [t]'
        else if( master == name )
            result << ' (default)'

        return result.toString()
    }

    @Deprecated
    List<String> getUpdates(int level) {

        def remote = strategy.lsRemote(false)
        List<String> branches = getBranchList()
                .findAll { it.name.startsWith('refs/heads/') || isRemoteBranch(it) }
                .unique { shortenRefName(it.name) }
                .findAll { Ref ref -> hasRemoteChange(ref,remote) }
                .collect { Ref ref -> formatUpdate(remote.get(ref.name),level) }

        remote = strategy.lsRemote(true)
        List<String> tags = getTagList()
                .findAll  { it.name.startsWith('refs/tags/') }
                .findAll { Ref ref -> hasRemoteChange(ref,remote) }
                .collect { Ref ref -> formatUpdate(remote.get(ref.name),level) }

        def result = new ArrayList<String>(branches.size() + tags.size())
        result.addAll(branches)
        result.addAll(tags)
        return result
    }

    /**
     * Checkout a specific revision, and fetch remote if not locally available.
     * @param revision The revision to be checked out
     */
    void checkout(String revision = null) {
        try {
            tryCheckout(revision)
        }
        catch( RefNotFoundException e ) {
            checkoutRemoteBranch(revision)
        }

    }

    /**
     * Checkout a specific revision returns exception if revision is not found locally
     * @param revision The revision to be checked out
     */
    void tryCheckout(String revision = null) throws RefNotFoundException {
        assert project
        if( !strategy ) {
            throw new IllegalStateException("Strategy not initialized")
        }
        strategy.tryCheckout(revision, manifest)
    }


    protected Ref checkoutRemoteBranch(String revision) {
        assert project
        if( !strategy ) {
            throw new IllegalStateException("Strategy not initialized")
        }
        return strategy.checkoutRemoteBranch(revision, manifest)
    }

    void updateModules() {

        if( !localPath )
            return // nothing to do

        final marker = new File(localPath, '.gitmodules')
        if( !marker.exists() || marker.empty() )
            return

        // the `gitmodules` attribute in the manifest makes it possible to enable/disable modules updating
        final modules = getManifest().gitmodules
        if( modules == false )
            return

        List<String> filter = []
        if( modules instanceof List<String> ) {
            filter.addAll(modules)
        }
        else if( modules instanceof String ) {
            filter.addAll( modules.tokenize(', ') )
        }

        final init = git.submoduleInit()
        final update = git.submoduleUpdate()

        if( manifest.recurseSubmodules ) {
            update.setStrategy(MergeStrategy.RECURSIVE)
        }

        filter.each { String m -> init.addPath(m); update.addPath(m) }
        // call submodule init
        init.call()
        // call submodule update
        if( provider.hasCredentials() )
            update.setCredentialsProvider( provider.getGitCredentials() )
        def updatedList = update.call()
        log.debug "Update submodules $updatedList"
    }

    protected String getRemoteCommitId(RevisionInfo rev) {
        final tag = rev.type == RevisionInfo.Type.TAG
        final list = strategy.lsRemote(tag).values()
        final ref = list.find { Repository.shortenRefName(it.name) == rev.name }
        if( !ref ) {
            log.debug "WARN: Cannot find any Git revision matching: ${rev.name}; ls-remote: $list"
            return null
        }
        return ref.objectId.name
    }

    protected void checkRemoteStatus0(RevisionInfo rev) {

        final remoteObjectId = getRemoteCommitId(rev)

        if( !remoteObjectId || remoteObjectId == rev.commitId ) {
            // nothing to do
            return
        }

        def local = rev.commitId.substring(0,10)
        def remote = remoteObjectId.substring(0,10)
        if( local == remote ) {
            remote = remoteObjectId
        }

        log.info "NOTE: Your local project version looks outdated - a different revision is available in the remote repository [$remote]"
    }

    void checkRemoteStatus(RevisionInfo rev) {
        try {
            checkRemoteStatus0(rev)
        }
        catch( Exception e ) {
            log.debug "WARN: Failed to check remote Git revision", e
        }
    }

    protected String getGitConfigRemoteUrl() {
        if( !isLocal() ) {
            return null
        }

        final gitConfig = localGitConfig
        if( !gitConfig.exists() ) {
            return null
        }

        final iniFile = new IniFile().load(gitConfig)
        final branch = manifest.getDefaultBranch()
        final remote = iniFile.getString("branch \"${branch}\"", "remote", "origin")
        final url = iniFile.getString("remote \"${remote}\"", "url")
        log.debug "Git config: $gitConfig; branch: $branch; remote: $remote; url: $url"
        return url
    }

    protected String getGitConfigRemoteDomain() {

        def str = getGitConfigRemoteUrl()
        if( !str ) return null

        try {
            final url = new GitUrl(str)
            if( url.protocol == 'file' ) {
                final hub = "file:${url.domain}"
                providerConfigs << new ProviderConfig(hub, [path:url.domain])
            }
            return url.domain
        }
        catch( IllegalArgumentException e) {
            log.debug e.message ?: e.toString()
            return null
        }

    }

    /**
     * Drop local copy of a repository. If revision is specified, only removes the specified revision
     * @param revision
     */
    void drop(String revision = null, boolean force = false) {
        if( isNotInitialized() )
            throw new AbortOperationException("No match found for: ${getProjectWithRevision()}")
        strategy.drop(revision, force)
    }


    protected String guessHubProviderFromGitConfig(boolean failFast=false) {
        assert localPath

        // find the repository remote URL from the git project config file
        final domain = getGitConfigRemoteDomain()
        if( !domain && failFast ) {
            def message = (localGitConfig.exists()
                            ? "Can't find git repository remote host -- Check config file at path: $localGitConfig"
                            : "Can't find git repository config file -- Repository may be corrupted: $localPath" )
            throw new AbortOperationException(message)
        }

        final result = domain ? providerConfigs.find { it -> it.domain == domain } : (ProviderConfig)null
        if( !result && failFast ) {
            def message = "Can't find any configured provider for git server `$domain` -- Make sure to have specified it in your `scm` file. For details check https://www.nextflow.io/docs/latest/sharing.html#scm-configuration-file"
            throw new AbortOperationException(message)
        }

        return result ? result.name : null
    }

    /**
     * Models revision information of a project repository.
     */
    @ToString
    @EqualsAndHashCode
    @TupleConstructor
    static class RevisionInfo {
        static enum Type {
            TAG,
            COMMIT,
            BRANCH
        }

        /**
         * Git commit ID
         */
        String commitId

        /**
         * Git tag or branch name
         */
        String name

        /**
         * The revision type.
         */
        Type type

        /**
         * @return A formatted string containing the commitId and revision properties
         */
        String toString() {

            if( !commitId ) {
                return '(unknown)'
            }

            if( name ) {
                return "${commitId.substring(0,10)} [${name}]"
            }

            return commitId
        }
    }

    /**
     * Enumeration of possible repository strategy types.
     */
    enum RepositoryStrategyType {
        /**
         * Legacy full clone strategy
         */
        LEGACY,

        /**
         * Bare repo + Multi-revision strategy
         */
        MULTI_REVISION
    }

}
