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
import static nextflow.util.StringUtils.redactUrlPassword

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import nextflow.config.ConfigParserFactory
import nextflow.config.Manifest
import nextflow.exception.AbortOperationException
import nextflow.exception.AmbiguousPipelineNameException
import nextflow.scm.HubOptions
import nextflow.script.ScriptFile
import nextflow.util.HttpClientOpts
import nextflow.util.IniFile
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.api.errors.RefNotFoundException
import org.eclipse.jgit.lib.Ref
import org.eclipse.jgit.lib.Repository
import org.eclipse.jgit.merge.MergeStrategy
/**
 * Handles operation on remote and local installed pipelines.
 * Local installations are managed by {@link MultiRevisionRepository}, which allows multiple
 * revisions to coexist efficiently by sharing objects through a bare repository and creating
 * lightweight clones for each commit.
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

    private String mainScript

    private RepositoryProvider provider

    private String hub

    private List<ProviderConfig> providerConfigs

    /**
     * The local pipeline repository
     */
    private MultiRevisionRepository repository

    /**
     * The http client settings declared in the SCM config file
     */
    private HttpClientOpts httpClientOpts

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
    AssetManager( String pipelineName, HubOptions hubOpts = null) {
        assert pipelineName
        // read the default config file (if available)
        def config = ProviderConfig.getDefault()
        // build the object
        build(pipelineName, config, hubOpts)
    }

    AssetManager( String pipelineName, Map config ) {
        assert pipelineName
        // build the object
        build(pipelineName, config)
    }

    AssetManager( String pipelineName, String revision, String mainScript = null, HubOptions hubOpts = null ) {
        assert pipelineName
        // read the default config file (if available)
        def config = ProviderConfig.getDefault()
        // build the object
        build(pipelineName, config, hubOpts, revision, mainScript)
    }

    /**
     * Build the asset manager internal data structure
     *
     * @param pipelineName A project name or a project repository Git URL
     * @param config A {@link Map} holding the configuration properties defined in the {@link ProviderConfig#DEFAULT_SCM_FILE} file
     * @param hubOpts The git provider credentials. See {@link HubOptions}
     * @return The {@link AssetManager} object itself
     */
    @PackageScope
    AssetManager build( String pipelineName, Map config = null, HubOptions hubOpts = null, String revision = null, String mainScript = null ) {

        this.providerConfigs = ProviderConfig.createFromMap(config)
        this.httpClientOpts = RepositoryProvider.httpClientOpts(config)

        this.project = resolveName(pipelineName)
        if( mainScript )
            this.mainScript = mainScript

        if( !isValidProjectName(project) )
            throw new IllegalArgumentException("Not a valid project name: ${project}")

        initRepository(revision)
        checkLegacyClone()
        this.hub = checkHubProvider(hubOpts)
        this.provider = createHubProvider(hub)

        if( revision )
            setRevision(revision)

        repository.setProvider(this.provider)

        setupCredentials(hubOpts)

        validateProjectDir()

        return this
    }

    private void initRepository(String revision = null) {
        repository = new MultiRevisionRepository(project, revision)
        if( revision )
            setRevision(revision)
        repository.setProvider(provider)
    }

    /**
     * Warn about a leftover clone created by the legacy asset layout, which is no longer used
     */
    private void checkLegacyClone() {
        final legacyPath = new File(root, project)
        if( new File(legacyPath, '.git').exists() )
            log.warn1 "The legacy pipeline cache at ${legacyPath} is no longer used and can be deleted (or migrated beforehand with `nextflow pull -migrate` using Nextflow 26.04 or 26.10)"
    }

    /**
     * Check if the local repository exists.
     * @return 'true' if it doesn't exist, otherwise 'false'.
     */
    boolean isNotInitialized() {
        return !MultiRevisionRepository.checkProject(root, project)
    }

    /**
     * Set the revision to work with
     */
    AssetManager setRevision(String revision) {
        if( provider )
            provider.setRevision(revision)
        if( revision && repository )
            repository.setRevision(revision)
        return this
    }

    /**
     * Get the current revision
     */
    private String getRevision() {
        if( provider )
            return provider.getRevision()
        return repository?.getRevision()
    }

    /**
     * Get the project name with revision appended
     */
    String getProjectWithRevision() {
        def rev = getRevision()
        return project + (rev ? ':' + rev : '')
    }

    @PackageScope
    File getLocalGitConfig() {
        return repository?.getLocalGitConfig()
    }

    @PackageScope
    AssetManager setProject(String name) {
        this.project = name
        if( !repository )
            initRepository(getRevision())
        repository.setProject(name)
        return this
    }

    protected RepositoryProvider getProvider() {
        return provider
    }

    /**
     * Sets the user credentials on the {@link RepositoryProvider} object
     *
     * @param hubOpts The git provider credentials. See {@link HubOptions}
     */
    @PackageScope
    void setupCredentials( HubOptions hubOpts ) {
        if( hubOpts?.getUser() ) {
            // rebind to the resolved hub provider so the password prompt names it correctly
            hubOpts = hubOpts.withProvider(hub)
            final user = hubOpts.getUser()
            final pwd = hubOpts.getPassword()
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

        if( !localPath?.exists() ) {
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
     * @param hubOpts The git provider credentials. See {@link HubOptions}
     * @return The name of hub name e.g. {@code github}, {@code bitbucket}, etc.
     */
    @PackageScope
    String checkHubProvider( HubOptions hubOpts ) {

        def result = hub
        if( !result )
            result = hubOpts?.getProvider()
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
        if( last.endsWith('.nf') ) {
            if( parts.size()==1 )
                throw new AbortOperationException("Not a valid project name: $name")

            if( mainScript )
                throw new AbortOperationException("Not a valid project name: $name -- Project name must be a directory when main script is provided")

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

        if( repository.endsWith('.nf') )
            throw new AbortOperationException("Repository URL must not end with a script file extension (.nf) -- use `-main-script` to specify the relative script path")

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

        return RepositoryFactory
            .newRepositoryProvider(config, project)
            .setHttpClientOpts(httpClientOpts ?: RepositoryProvider.httpClientOpts())
    }

    AssetManager setLocalPath(File path) {
        if( !repository )
            initRepository(getRevision())
        repository.setLocalPath(path)
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
        return repository?.getGitRepositoryUrl() ?: provider.getCloneUrl()
    }

    File getLocalPath() {
        return repository?.getLocalPath()
    }

    File getProjectPath() {
        return repository?.getProjectPath()
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
        if( !localPath?.exists() ) {
            throw new AbortOperationException("Unknown project folder: ${localPath ?: project}")
        }

        def mainScript = scriptName ?: getMainScriptName()
        def result = new File(localPath, mainScript)
        if( !result.exists() )
            throw new AbortOperationException("Missing project main script: ${result.canonicalPath}")

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
        return repository?.getRemoteDefaultBranch() ?: DEFAULT_BRANCH
    }

    @Memoized
    Manifest getManifest() {
        getManifest0()
    }


    protected Manifest getManifest0() {
        String text = null
        ConfigObject result = null
        try {
            text = localPath?.exists() ? new File(localPath, MANIFEST_FILE_NAME).text : provider.readText(MANIFEST_FILE_NAME)
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
        if( localPath?.exists() ) {
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
     * @return {@code true} when the SCM source points to a local file-system repository.
     */
    boolean isLocalScmSource() {
        return provider instanceof LocalRepositoryProvider
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
        return repository?.isClean() ?: true
    }

    /**
     * Close the underlying Git repository
     */
    void close() {
        repository?.close()
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
     * Append the projects found in the '.repos' subdirectory
     * @param dir
     * @param result
     */
    private static void appendProjectsFromDir(File dir, HashSet result) {
        final reposDir = new File(dir, MultiRevisionRepository.REPOS_SUBDIR)
        if( !reposDir.isDirectory() )
            return
        reposDir.eachDir { File org ->
            org.eachDir { File it ->
                log.debug("Adding ${org.getName()}/${it.getName()}")
                result << "${org.getName()}/${it.getName()}".toString()
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
        return repository?.getGit()
    }

    /**
     * Download a pipeline from a remote Github repository
     *
     * @param revision The revision to download
     * @return A message representing the operation result
     */
    String download(String revision = null) {
        assert project
        if( !repository ) {
            throw new IllegalStateException("Repository not initialized")
        }
        // If it is a new download check is a valid repository
        if( !localPath?.exists() ) {
            checkValidRemoteRepo(revision)
        }
        return repository.download(revision, manifest)
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
        log.debug "Clone project `$project` -- Using remote URI: ${redactUrlPassword(uri)} into: $directory"

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
        if( !repository ) {
            throw new IllegalStateException("Repository not initialized")
        }
        return repository.getCurrentRevision()
    }

    RevisionInfo getCurrentRevisionAndName() {
        if( !repository ) {
            throw new IllegalStateException("Repository not initialized")
        }
        return repository.getCurrentRevisionAndName()
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

        def pulled = repository.listDownloadedCommits()
        def master = getDefaultBranch()

        List<String> branches = getBranchList()
            .findAll { it.name.startsWith('refs/heads/') || isRemoteBranch(it) }
            .unique { shortenRefName(it.name) }
            .collect { Ref it -> refToString(it, pulled, master, false, level) }

        List<String> tags = getTagList()
                .findAll  { it.name.startsWith('refs/tags/') }
                .collect { refToString(repository.peel(it) , pulled, master, true, level) }

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
        final currentCommits = repository.listDownloadedCommits()
        final master = getDefaultBranch()

        final branches = []
        final tags = []
        final pulled = new LinkedList<Map>()
        log.debug("Current commits : $currentCommits")
        Map<String, Ref> remote = checkForUpdates ? repository.lsRemote(false) : null
        getBranchList()
            .findAll { it.name.startsWith('refs/heads/') || isRemoteBranch(it) }
            .unique { shortenRefName(it.name) }
            .each {
                final map = refToMap(it, remote)
                if( isRefInCommits(it, currentCommits) )
                    pulled << map
                branches << map
            }

        remote = checkForUpdates ? repository.lsRemote(true) : null
        getTagList()
            .findAll { it.name.startsWith('refs/tags/') }
            .each {
                Ref ref = repository.peel(it)
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
        return repository?.getBranchList() ?: []
    }

    @Memoized
    protected List<Ref> getTagList() {
        return repository?.getTagList() ?: []
    }

    protected String refToString(Ref ref, List<String> downloaded, String master, boolean tag, int level ) {

        def result = new StringBuilder()
        def name = shortenRefName(ref.name)
        result << ( isRefInCommits(ref, downloaded) ? '>' : ' ')

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

        def remote = repository.lsRemote(false)
        List<String> branches = getBranchList()
                .findAll { it.name.startsWith('refs/heads/') || isRemoteBranch(it) }
                .unique { shortenRefName(it.name) }
                .findAll { Ref ref -> hasRemoteChange(ref,remote) }
                .collect { Ref ref -> formatUpdate(remote.get(ref.name),level) }

        remote = repository.lsRemote(true)
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
        if( !repository ) {
            throw new IllegalStateException("Repository not initialized")
        }
        repository.tryCheckout(revision)
    }


    protected Ref checkoutRemoteBranch(String revision) {
        assert project
        if( !repository ) {
            throw new IllegalStateException("Repository not initialized")
        }
        return repository.checkoutRemoteBranch(revision, manifest)
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
        if( modules instanceof List ) {
            filter.addAll(modules as List<String>)
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
        final list = repository.lsRemote(tag).values()
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
        if( !gitConfig?.exists() ) {
            return null
        }

        final iniFile = new IniFile().load(gitConfig)
        final url = iniFile.getString('remote "origin"', "url")
        log.debug "Git config: $gitConfig; url: ${redactUrlPassword(url)}"
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
        repository.drop(revision, force)
    }


    protected String guessHubProviderFromGitConfig(boolean failFast=false) {
        // find the repository remote URL from the git project config file
        final domain = getGitConfigRemoteDomain()
        if( !domain && failFast ) {
            def message = (localGitConfig?.exists()
                            ? "Can't find git repository remote host -- Check config file at path: $localGitConfig"
                            : "Can't find git repository config file -- Repository may be corrupted: $localPath" )
            throw new AbortOperationException(message)
        }

        final result = domain ? providerConfigs.find { it -> it.domain == domain } : (ProviderConfig)null
        if( !result && failFast ) {
            def message = "Can't find any configured provider for git server `$domain` -- Make sure to have specified it in your `scm` file. For details check https://docs.seqera.io/nextflow/sharing#git-configuration"
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

}
