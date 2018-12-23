/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import nextflow.cli.HubOptions
import nextflow.config.ConfigParser
import nextflow.config.Manifest
import nextflow.exception.AbortOperationException
import nextflow.script.ScriptFile
import nextflow.util.IniFile
import org.eclipse.jgit.api.CreateBranchCommand
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.api.ListBranchCommand
import org.eclipse.jgit.api.MergeResult
import org.eclipse.jgit.api.errors.RefNotFoundException
import org.eclipse.jgit.errors.RepositoryNotFoundException
import org.eclipse.jgit.lib.Constants
import org.eclipse.jgit.lib.ObjectId
import org.eclipse.jgit.lib.Ref
import org.eclipse.jgit.lib.Repository
import org.eclipse.jgit.transport.UsernamePasswordCredentialsProvider
import static nextflow.Const.DEFAULT_HUB
import static nextflow.Const.DEFAULT_MAIN_FILE_NAME
import static nextflow.Const.DEFAULT_ORGANIZATION
import static nextflow.Const.DEFAULT_ROOT
import static nextflow.Const.MANIFEST_FILE_NAME
/**
 * Handles operation on remote and local installed pipelines
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
class AssetManager {

    /**
     * The folder all pipelines scripts are installed
     */
    @PackageScope
    static File root = DEFAULT_ROOT

    /**
     * The pipeline name. It must be in the form {@code username/repo} where 'username'
     * is a valid user name or organisation account, while 'repo' is the repository name
     * containing the pipeline code
     */
    private String project

    /**
     * Directory where the pipeline is cloned (i.e. downloaded)
     */
    private File localPath

    private Git _git

    private String mainScript

    private RepositoryProvider provider

    private String hub

    private List<ProviderConfig> providerConfigs

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

    /**
     * Build the asset manager internal data structure
     *
     * @param pipelineName A project name or a project repository Git URL
     * @param config A {@link Map} holding the configuration properties defined in the {@link ProviderConfig#SCM_FILE} file
     * @param cliOpts User credentials provided on the command line. See {@link HubOptions} trait
     * @return The {@link AssetManager} object itself
     */
    @PackageScope
    AssetManager build( String pipelineName, Map config = null, HubOptions cliOpts = null ) {

        this.providerConfigs = ProviderConfig.createFromMap(config)

        this.project = resolveName(pipelineName)
        this.localPath = checkProjectDir(project)
        this.hub = checkHubProvider(cliOpts)
        this.provider = createHubProvider(hub)
        setupCredentials(cliOpts)
        validateProjectDir()

        return this
    }

    @PackageScope
    File getLocalGitConfig() {
        localPath ? new File(localPath,'.git/config') : null
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
     * Verify the project name matcher the expected pattern.
     * and return the directory where the project is stored locally
     *
     * @param projectName A project name matching the pattern {@code owner/project}
     * @return The project dir {@link File}
     */
    @PackageScope
    File checkProjectDir(String projectName) {

        if( !isValidProjectName(projectName)) {
            throw new IllegalArgumentException("Not a valid project name: $projectName")
        }

        new File(root, project)
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
            throw new AbortOperationException("A project with name: `$localPath` has been already download from a different provider: `$configProvider`")
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
    String resolveName( String name ) {
        assert name

        def project = checkForGitUrl(name)
        if( project )
            return project

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
            throw new AbortOperationException("Which one do you mean?\n${qualifiedName.join('\n')}")
        }

        return qualifiedName
    }

    String getProject() { project }

    String getHub() { hub }

    @PackageScope
    String checkForGitUrl( String repository ) {

        if( repository.startsWith('http://') || repository.startsWith('https://') || repository.startsWith('file:/')) {
            try {
                def url = new GitUrl(repository)

                def result
                if( url.protocol == 'file' ) {
                    this.hub = "file:${url.domain}"
                    providerConfigs << new ProviderConfig(this.hub, [path:url.domain])
                    result = "local/${url.project}"
                }
                else {
                    // find the provider config for this server
                    this.hub = providerConfigs.find { it.domain == url.domain } ?.name
                    result = url.project
                }
                log.debug "Repository URL: $repository; Project: $result; Hub provider: $hub"

                return result
            }
            catch( IllegalArgumentException e ) {
                log.debug "Cannot parse Git URL: $repository -- cause: ${e.message}"
            }
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

        RepositoryProvider .create(config, project)

    }

    AssetManager setLocalPath(File path) {
        this.localPath = path
        return this
    }

    AssetManager setForce( boolean value ) {
        this.force = value
        return this
    }

    AssetManager checkValidRemoteRepo() {
        def scriptName = getMainScriptName()
        provider.validateFor(scriptName)
        return this
    }

    @Memoized
    String getGitRepositoryUrl() {

        if( localPath.exists() ) {
            return localPath.toURI().toString()
        }

        provider.getCloneUrl()
    }

    File getLocalPath() { localPath }

    ScriptFile getScriptFile() {

        def result = new ScriptFile(getMainScriptFile())
        result.revisionInfo = getCurrentRevisionAndName()
        result.repository = getGitConfigRemoteUrl()
        result.localPath = localPath.toPath()

        return result
    }

    File getMainScriptFile() {
        if( !localPath.exists() ) {
            throw new AbortOperationException("Unknown project folder: $localPath")
        }

        def mainScript = getMainScriptName()
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
        getManifest().getDefaultBranch()
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
            def config = new ConfigParser().setIgnoreIncludes(true).parse(text)
            result = (ConfigObject)config.manifest
        }
        catch( Exception e ) {
            throw new Exception("Project config file is malformed -- Cause: ${e.message ?: e}", e)
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
            catch (IOException e) {
                provider.validateRepo()
                log.debug "Cannot retried remote config file -- likely does not exist"
                return null
            }
        }
    }


    String getBaseName() {
        def result = project.tokenize('/')
        if( result.size() > 2 ) throw new IllegalArgumentException("Not a valid projct name: $project")
        return result.size()==1 ? result[0] : result[1]
    }

    boolean isLocal() {
        localPath.exists()
    }

    /**
     * @return {@code true} when the local project path exists and contains at least the default script
     *      file (i.e. main.nf) or the nextflow manifest file (i.e. nextflow.config)
     */
    boolean isRunnable() {
        localPath.exists() && ( new File(localPath,DEFAULT_MAIN_FILE_NAME).exists() || new File(localPath,MANIFEST_FILE_NAME).exists() )
    }

    /**
     * @return true if no differences exist between the working-tree, the index,
     *         and the current HEAD, false if differences do exist
     */
    boolean isClean() {
        try {
            git.status().call().isClean()
        }
        catch( RepositoryNotFoundException e ) {
            return true
        }
    }

    /**
     * Close the underlying Git repository
     */
    void close() {
        if( _git ) {
            _git.close()
            _git = null
        }
    }

    /**
     * @return The list of available pipelines
     */
    static List<String> list() {
        log.debug "Listing projects in folder: $root"

        def result = []
        if( !root.exists() )
            return result

        root.eachDir { File org ->
            org.eachDir { File it ->
                result << "${org.getName()}/${it.getName()}".toString()
            }
        }

        return result
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
        if( !_git ) {
            _git = Git.open(localPath)
        }
        return _git
    }

    /**
     * Download a pipeline from a remote Github repository
     *
     * @param revision The revision to download
     * @result A message representing the operation result
     */
    def download(String revision=null) {
        assert project

        /*
         * if the pipeline already exists locally pull it from the remote repo
         */
        if( !localPath.exists() ) {
            localPath.parentFile.mkdirs()
            // make sure it contains a valid repository
            checkValidRemoteRepo()

            final cloneURL = getGitRepositoryUrl()
            log.debug "Pulling $project -- Using remote clone url: ${cloneURL}"

            // clone it
            def clone = Git.cloneRepository()
            if( provider.hasCredentials() )
                clone.setCredentialsProvider( new UsernamePasswordCredentialsProvider(provider.user, provider.password) )

            if( revision ) {
                clone.setBranch(revision)
            }

            clone
                .setURI(cloneURL)
                .setDirectory(localPath)
                .call()

            // return status message
            return "downloaded from ${cloneURL}"
        }


        log.debug "Pull pipeline $project  -- Using local path: $localPath"

        // verify that is clean
        if( !isClean() )
            throw new AbortOperationException("$project contains uncommitted changes -- cannot pull from repository")

        if( revision && revision != getCurrentRevision() ) {
            /*
             * check out a revision before the pull operation
             */
            try {
                git.checkout() .setName(revision) .call()
            }
            /*
             * If the specified revision does not exist
             * Try to checkout it from a remote branch and return
             */
            catch ( RefNotFoundException e ) {
                def ref = checkoutRemoteBranch(revision)
                return "checkout-out at ${ref.getObjectId()}"
            }
        }

        def pull = git.pull()
        def revInfo = getCurrentRevisionAndName()

        if ( revInfo.type == RevisionInfo.Type.COMMIT ) {
            log.debug("Repo appears to be checked out to a commit hash, but not a TAG, so we are NOT pulling the repo and assuming it is already up to date!")
            return MergeResult.MergeStatus.ALREADY_UP_TO_DATE.toString()
        }

        if ( revInfo.type == RevisionInfo.Type.TAG ) {
            pull.setRemoteBranchName( "refs/tags/" + revInfo.name )
        }

        if( provider.hasCredentials() )
            pull.setCredentialsProvider( new UsernamePasswordCredentialsProvider(provider.user, provider.password))

        def result = pull.call()
        if(!result.isSuccessful())
            throw new AbortOperationException("Cannot pull project `$project` -- ${result.toString()}")

        return result?.mergeResult?.mergeStatus?.toString()

    }

    /**
     * Clone a pipeline from a remote pipeline repository to the specified folder
     *
     * @param directory The folder when the pipeline will be cloned
     * @param revision The revision to be cloned. It can be a branch, tag, or git revision number
     */
    void clone(File directory, String revision = null) {

        def clone = Git.cloneRepository()
        def uri = getGitRepositoryUrl()
        log.debug "Clone project `$project` -- Using remote URI: ${uri} into: $directory"

        if( !uri )
            throw new AbortOperationException("Cannot find the specified project: $project")

        clone.setURI(uri)
        clone.setDirectory(directory)
        if( provider.hasCredentials() )
            clone.setCredentialsProvider(new UsernamePasswordCredentialsProvider(provider.user, provider.password))

        if( revision )
            clone.setBranch(revision)

        clone.call()
    }

    /**
     * @return The symbolic name of the current revision i.e. the current checked out branch or tag
     */
    String getCurrentRevision() {
        Ref head = git.getRepository().findRef(Constants.HEAD);
        if( !head )
            return '(unknown)'

        if( head.isSymbolic() )
            return Repository.shortenRefName(head.getTarget().getName())

        if( !head.getObjectId() )
            return '(unknown)'

        // try to resolve the the current object id to a tag name
        Map<ObjectId, String> names = git.nameRev().addPrefix( "refs/tags/" ).add(head.objectId).call()
        names.get( head.objectId ) ?: head.objectId.name()
    }

    RevisionInfo getCurrentRevisionAndName() {
        Ref head = git.getRepository().findRef(Constants.HEAD);
        if( !head )
            return null

        if( head.isSymbolic() ) {
            return new RevisionInfo(head.objectId.name(), Repository.shortenRefName(head.getTarget().getName()), RevisionInfo.Type.BRANCH)
        }

        if( !head.getObjectId() )
            return null

        // try to resolve the the current object id to a tag name
        Map<ObjectId, String> allNames = git.nameRev().addPrefix( "refs/tags/" ).add(head.objectId).call()
        def name = allNames.get( head.objectId )
        if( name ) {
            return new RevisionInfo(head.objectId.name(), name, RevisionInfo.Type.TAG)
        }
        else {
            return new RevisionInfo(head.objectId.name(), null, RevisionInfo.Type.COMMIT)
        }
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

        def current = getCurrentRevision()
        def master = getDefaultBranch()

        List<String> branches = getBranchList()
            .findAll { it.name.startsWith('refs/heads/') || it.name.startsWith('refs/remotes/origin/') }
            .unique { shortenRefName(it.name) }
            .collect { Ref it -> refToString(it,current,master,false,level) }

        List<String> tags = getTagList()
                .findAll  { it.name.startsWith('refs/tags/') }
                .collect { refToString(it,current,master,true,level) }

        def result = new ArrayList<String>(branches.size() + tags.size())
        result.addAll(branches)
        result.addAll(tags)
        return result
    }

    Map getBranchesAndTags(boolean checkForUpdates) {
        final result = [:]
        final current = getCurrentRevision()
        final master = getDefaultBranch()

        final branches = []
        final tags = []

        Map<String, Ref> remote = checkForUpdates ? git.lsRemote().callAsMap() : null
        getBranchList()
                .findAll { it.name.startsWith('refs/heads/') || it.name.startsWith('refs/remotes/origin/') }
                .unique { shortenRefName(it.name) }
                .each { Ref it -> branches << refToMap(it,remote)  }

        remote = checkForUpdates ? git.lsRemote().setTags(true).callAsMap() : null
        getTagList()
                .findAll  { it.name.startsWith('refs/tags/') }
                .each { Ref it -> tags << refToMap(it,remote) }

        result.current = current    // current branch name
        result.master = master      // master branch name
        result.branches = branches  // collection of branches
        result.tags = tags          // collect of tags 
        return result
    }

    protected Map refToMap(Ref ref, Map<String,Ref> remote) {
        final entry = new HashMap(2)
        final peel = git.getRepository().peel(ref)
        final objId = peel.getPeeledObjectId() ?: peel.getObjectId()
        // the branch or tag name
        entry.name = shortenRefName(ref.name)
        // the local commit it
        entry.commitId = objId.name()
        // the remote commit Id for this branch or tag
        if( remote && hasRemoteChange(ref,remote) ) {
            entry.latestId = remote.get(ref.name).objectId.name()
        }
        return entry
    }

    @Memoized
    protected List<Ref> getBranchList() {
        git.branchList().setListMode(ListBranchCommand.ListMode.ALL) .call()
    }

    @Memoized
    protected List<Ref> getTagList() {
        git.tagList().call()
    }

    protected formatObjectId(ObjectId obj, boolean human) {
        return human ? obj.name.substring(0,10) : obj.name
    }

    protected String refToString(Ref ref, String current, String master, boolean tag, int level ) {

        def result = new StringBuilder()
        def name = shortenRefName(ref.name)
        result << (name == current ? '*' : ' ')

        if( level ) {
            def peel = git.getRepository().peel(ref)
            def obj = peel.getPeeledObjectId() ?: peel.getObjectId()
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

    private String shortenRefName( String name ) {
        if( name.startsWith('refs/remotes/origin/') )
            return name.replace('refs/remotes/origin/', '')

        return Repository.shortenRefName(name)
    }

    protected String formatUpdate(Ref remoteRef, int level) {

        def result = new StringBuilder()
        result << '  '
        result << formatObjectId(remoteRef.objectId, level<2)
        result << ' '
        result << shortenRefName(remoteRef.name)

        return result.toString()
    }

    protected hasRemoteChange(Ref ref, Map<String,Ref> remote) {
        if( !remote.containsKey(ref.name) )
            return false
        ref.objectId.name != remote[ref.name].objectId.name
    }

    @Deprecated
    List<String> getUpdates(int level) {

        def remote = git.lsRemote().callAsMap()
        List<String> branches = getBranchList()
                .findAll { it.name.startsWith('refs/heads/') || it.name.startsWith('refs/remotes/origin/') }
                .unique { shortenRefName(it.name) }
                .findAll { Ref ref -> hasRemoteChange(ref,remote) }
                .collect { Ref ref -> formatUpdate(remote.get(ref.name),level) }

        remote = git.lsRemote().setTags(true).callAsMap()
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
     * Checkout a specific revision
     * @param revision The revision to be checked out
     */
    void checkout( String revision = null ) {
        assert localPath

        def current = getCurrentRevision()
        if( current != defaultBranch ) {
            if( !revision ) {
                throw new AbortOperationException("Project `$project` currently is sticked on revision: $current -- you need to specify explicitly a revision with the option `-r` to use it")
            }
        }
        else if( !revision || revision == current ) {
            // nothing to do
            return
        }

        // verify that is clean
        if( !isClean() )
            throw new AbortOperationException("Project `$project` contains uncommitted changes -- Cannot switch to revision: $revision")

        try {
            git.checkout().setName(revision) .call()
        }
        catch( RefNotFoundException e ) {
            checkoutRemoteBranch(revision)
        }

    }


    protected Ref checkoutRemoteBranch( String revision ) {

        try {
            git.fetch().call()
            git.checkout()
                    .setCreateBranch(true)
                    .setName(revision)
                    .setUpstreamMode(CreateBranchCommand.SetupUpstreamMode.TRACK)
                    .setStartPoint("origin/" + revision)
                    .call()
        }
        catch (RefNotFoundException e) {
            throw new AbortOperationException("Cannot find revision `$revision` -- Make sure that it exists in the remote repository `$repositoryUrl`", e)
        }

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
            filter.addAll(modules as List)
        }
        else if( modules instanceof String ) {
            filter.addAll( (modules as String).tokenize(', ') )
        }

        final init = git.submoduleInit()
        final update = git.submoduleUpdate()
        filter.each { String m -> init.addPath(m); update.addPath(m) }
        // call submodule init
        init.call()
        // call submodule update
        def updatedList = update.call()
        log.debug "Update submodules $updatedList"
    }

    protected String getRemoteCommitId(RevisionInfo rev) {
        final tag = rev.type == RevisionInfo.Type.TAG
        final list = git.lsRemote().setTags(tag).call()
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
        if( !localPath ) {
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


    protected String guessHubProviderFromGitConfig(boolean failFast=false) {
        assert localPath
        
        // find the repository remote URL from the git project config file
        final server = getGitConfigRemoteDomain()
        if( !server && failFast ) {
            def message = (localGitConfig.exists()
                            ? "Can't find git repository remote host -- Check config file at path: $localGitConfig"
                            : "Can't find git repository config file -- Repository may be corrupted: $localPath" )
            throw new AbortOperationException(message)
        }

        final result = providerConfigs.find { it -> it.domain == server }
        if( !result && failFast ) {
            def message = "Can't find any configured provider for git server `$server` -- Make sure to have specified it in your `scm` file. For details check https://www.nextflow.io/docs/latest/sharing.html#scm-configuration-file"
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

            commitId
        }
    }
}
