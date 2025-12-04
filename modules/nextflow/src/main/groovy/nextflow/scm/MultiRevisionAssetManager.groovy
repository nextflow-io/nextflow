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
 */

package nextflow.scm

import nextflow.exception.AbortOperationException
import nextflow.file.FileMutex
import org.eclipse.jgit.internal.storage.file.FileRepository
import org.eclipse.jgit.storage.file.FileRepositoryBuilder
import org.eclipse.jgit.transport.RefSpec

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.cli.HubOptions
import nextflow.util.IniFile
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.api.ListBranchCommand
import org.eclipse.jgit.errors.RepositoryNotFoundException
import org.eclipse.jgit.lib.Constants
import org.eclipse.jgit.lib.Ref
import org.eclipse.jgit.lib.Repository
/**
 * Handles operation on remote and local installed pipelines
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
class MultiRevisionAssetManager extends AssetManager{

    static final String BARE_REPO = '.nextflow/bare_repo'

    static public final String REVISION_SUBDIR = '.nextflow/commits'

    /**
     * The name of the commit/branch/tag as requested via command line
     * This is now a first class attribute of a pipeline
     */
    private String revision

    private File commitPath

    private Git _bareGit

    private Git _commitGit

    /**
     * Create a new asset manager object with default parameters
     */
    MultiRevisionAssetManager() {
        super()
    }

    /**
     * Create a new asset manager with the specified pipeline name
     *
     * @param pipelineName The pipeline to be managed by this manager e.g. {@code nextflow-io/hello}
     */
    MultiRevisionAssetManager(String pipelineName, HubOptions cliOpts = null, String revision= null) {
        assert pipelineName
        // read the default config file (if available)
        def config = ProviderConfig.getDefault()
        // build the object
        build(pipelineName, config, cliOpts, revision)
    }

    MultiRevisionAssetManager(String pipelineName, Map config, String revision= null) {
        assert pipelineName
        // build the object
        build(pipelineName, config)
    }

    @PackageScope
    MultiRevisionAssetManager build(String pipelineName, Map config = null, HubOptions cliOpts = null, String revision = null ) {
        super.build(pipelineName, config, cliOpts)
        setRevision(revision)
        return this
    }

    @Override
    File getLocalGitConfig() {
        bareRepo.exists() ? new File(bareRepo,'config')
            : super.getLocalGitConfig()
    }

    @PackageScope
    File getBareRepo() {
        new File(root, project + '/' + BARE_REPO)
    }

    @PackageScope
    File getRevisionSubdir( String projectName = project ) {
        new File(root, projectName + '/' + REVISION_SUBDIR)
    }

    /**
     * Update the directory where the commitId shared clone will be stored.
     *
     * @param revision Revision ID for the selected pipeline (git branch, tag or commit SHA number)
     * @return The project dir {@link File}
     */
    @PackageScope
    void updateCommitDir(String commitId) {
        if( !commitId )
            return
        final oldCommitPath = this.commitPath
        if( oldCommitPath && _commitGit) {
            _commitGit.close()
            _commitGit = null
        }
        this.commitPath = new File( root, project + '/' + REVISION_SUBDIR + '/' + commitId )
    }

    @Override
    String getGitRepositoryUrl() {
        provider.getCloneUrl()
    }

    @PackageScope
    void checkBareRepo() {
        /*
         * if the bare repository of the pipeline does not exists locally pull it from the remote repo
         */
        if( !bareRepo.exists() ) {
            bareRepo.parentFile.mkdirs()

            final cloneURL = getGitRepositoryUrl()
            log.debug "Pulling bare repo for $project -- Using remote clone url: ${cloneURL}"

            def bare = Git.cloneRepository()
            if( provider.hasCredentials() )
                bare.setCredentialsProvider( provider.getGitCredentials() )

            bare
                .setBare( true )
                .setURI(cloneURL)
                .setGitDir(bareRepo)
                .call()
        }  else  {
            final updateRevision = revision ?: getDefaultBranch()
            log.debug "Fetching (updating) bare repo for $project [revision: $updateRevision]"
            bareGit.fetch().setRefSpecs(refSpecForName(updateRevision)).call()
        }
    }

    @PackageScope
    String revisionToCommitWithBareRepo(String revision) {
        String commitId = null

        if( bareRepo.exists() ) {
            def rev = Git.open(bareRepo)
                         .getRepository()
                         .resolve(revision ?: Constants.HEAD)
            if( rev )
                commitId = rev.getName()
        }

        return commitId
    }


    MultiRevisionAssetManager setRevision(String revision) {
        if(! revision )
            revision = getDefaultBranchFromBareRepo()
            this.revision = revision
        updateCommitDir(revisionToCommitWithBareRepo(revision))
        return this
    }

    String getRevision() { revision }

    String getProjectWithRevision() { project + ( revision ? ':' + revision : '' ) }

    @Override
    File getLocalPath() {
        return this.commitPath ?: super.getLocalPath()
    }

    @Override
    protected String getRemoteDefaultBranch() {
        return getDefaultBranchFromBareRepo()
            ?: super.getRemoteDefaultBranch()
    }

    protected String getDefaultBranchFromBareRepo(){
        if (!bareRepo.exists()) {
            log.debug "Bare repo ${bareRepo} doesn't exist"
            return null
        }
        def head = getDefaultBranchRef()
        if( head ) {
            Repository.shortenRefName(head.getName())
        } else {
            log.debug "Unable to determine default branch from bare repo ${bareRepo}"
            return null
        }
    }

    protected Ref getDefaultBranchRef(){
        return bareGit.getRepository().findRef(Constants.HEAD).getTarget()
    }

    /**
     * @return True if MultiRevisionAssetManager already contains the bare repository of the project
     */
    boolean hasBareRepo() {
        bareRepo.exists()
    }

    @Override
    boolean isClean() {
        try {
            commitGit.status().call().isClean()
        }
        catch( RepositoryNotFoundException e ) {
            return true
        }
    }

    /**
     * Close the underlying Git repository
     */
    void close() {
        super.close()
        if( _bareGit ) {
            _bareGit.close()
            _bareGit = null
        }
        if( _commitGit ) {
            _commitGit.close()
            _commitGit = null
        }
    }

    /**
     * @return The list of available revisions for a given project name
     */
    List<String> listRevisions( String projectName = this.project ) {
        log.debug "Listing revisions for project: $projectName"
        if( !root.exists() )
            return []
        return getBranchesAndTags(false).pulled as List<String> ?: []
    }

    /**
     * @return The list of downloaded bare commits for a given project name
     */
    List<String> listCommits( String projectName = this.project ) {
        log.debug "Listing all commits for project: $projectName"
        def result = new LinkedList()
        if( !root.exists() )
            return result
        if( !getRevisionSubdir(projectName).exists() )
            return result

        getRevisionSubdir(projectName).eachDir { File it -> result << it.getName().toString() }

        return result
    }

    protected Git getBareGit() {
        if( !_bareGit ) {
            _bareGit = Git.open(bareRepo)
        }
        return _bareGit
    }

    protected Git getCommitGit() {
        if( !_commitGit ) {
            _commitGit = Git.open(commitPath)
        }
        return _commitGit
    }

    @Override
    protected Git getGit(){
        if( commitPath && commitPath.exists() )
            return getCommitGit()
        return super.git
    }

    private Git getBareGitWithLegacyFallback() {
        if( hasBareRepo() )
            return bareGit
        return super.git
    }

    /**
     * Creates a clone in commits folder where objects are shared with the bare repository.
     *
     * @param testDisableUpdateLocalPath
     * @param createTimeout
     * @return
     */
    String createSharedClone(boolean testDisableUpdateLocalPath = false, createTimeout = '60s'){
        assert project

        // make sure it contains a valid repository
        checkValidRemoteRepo()

        // get local copy of bare repository
        checkBareRepo()
        // update mapping of revision to commit, and update localPath
        // boolean is for testing purposes only (e.g. UpdateModuleTest)
        final downloadRevision = revision ?: getDefaultBranch()
        final commitId = revisionToCommitWithBareRepo(downloadRevision)
        if (!commitId){
            throw new AbortOperationException("No commit found for revision ${downloadRevision} in project $project")
        }
        updateCommitDir(commitId)
        if( commitPath.exists() )
            return
        /*
         * if revision does not exists locally pull it from the remote repo
         */
        if( !commitPath.parentFile.exists() )
            commitPath.parentFile.mkdirs()

        // Use a file mutex to prevent concurrent clones of the same commit.
        final file = new File(commitPath.parentFile, ".${commitPath.name}.lock")
        final wait = "Another Nextflow instance is creating clone for $downloadRevision -- please wait till it completes"
        final err = "Unable to acquire exclusive lock after $createTimeout on file: $file"

        final mutex = new FileMutex(target: file, timeout: createTimeout, waitMessage: wait, errorMessage: err)
        try {
            mutex.lock { createSharedClone0(downloadRevision) }
        }
        finally {
            file.delete()
        }
        return "downloaded from ${getGitRepositoryUrl()}"
    }

    private String createSharedClone0(String downloadRevision) {
        // This check is required in case of two nextflow instances were doing the same shared clone at the same time.
        // If commitPath exists the previous nextflow instance created the shared clone succesfully.
        if( commitPath.exists() )
            return

        commitPath.parentFile.mkdirs()
        try {
            final cloneURL = bareRepo.toString()
            log.debug "Pulling $project -- Using remote clone url: ${cloneURL}"

            // clone it, but don't specify a revision - jgit will checkout the default branch
            File bareObjectsDir = new File(bareRepo, "objects")
            FileRepository repo = new FileRepositoryBuilder().setGitDir(new File(localPath, ".git")).addAlternateObjectDirectory(bareObjectsDir).build() as FileRepository
            repo.create()

            // Write alternates file (not done by repo create).
            new File(repo.getObjectsDirectory(), "info/alternates").write(bareObjectsDir.absolutePath)

            // 3️⃣ Configure remote pointing to the cache repo
            repo.getConfig().setString("remote", "origin", "url", cloneURL);
            repo.getConfig().save();


            commitGit.fetch().setRefSpecs(refSpecForName(downloadRevision)).call()

            commitGit.checkout().setName(downloadRevision).call()

        } catch (Throwable t){
            // If there is an error creating the sared clone, remove the local path to avoid incorrect clones.
            commitPath.deleteDir()
            throw t
        }

        // return status message
        return "downloaded from ${getGitRepositoryUrl()}"
    }

    private RefSpec refSpecForName(String revision) {
        // Is it a local branch?
        Ref branch = bareGit.getRepository().findRef("refs/heads/" + revision);
        if (branch != null) {
            return new RefSpec("refs/heads/" + revision + ":refs/heads/" + revision);
        }

        // Is it a tag?
        Ref tag = bareGit.getRepository().findRef("refs/tags/" + revision);
        if (tag != null) {
            return new RefSpec("refs/tags/" + revision + ":refs/tags/" + revision);
        }

        // It is a commit
        return new RefSpec(revision + ":refs/tags/" + revision)
    }

    /**
     * @return A list of existing branches and tags names. For example
     * <pre>
     *     * P master (default)
     *         patch-x
     *         v1.0 (t)
     *         v1.1 (t)
     * </pre>
     *
     * The character {@code P} on the left indicates the revision is pulled locally,
     *  the string {@code (default)} ticks that it is the default working branch,
     *  while the string {@code (t)} shows that the revision is a git tag (instead of a branch)
     */
    @Deprecated
    @Override
    List<String> getRevisions(int level) {

        def current = getCurrentRevision()
        def master = getDefaultBranch()
        def pulled = listCommits()
        def headRef = getDefaultBranchRef()
        List<String> branches = getBranchList()
            .findAll { (it.name.startsWith('refs/heads/') && it.name != headRef?.name ) || isRemoteBranch(it) }
            .unique { shortenRefName(it.name) }
            .collect { Ref it -> refToString(it,current,master,pulled,false,level) }

        List<String> tags = getTagList()
                .findAll  { it.name.startsWith('refs/tags/') }
                .collect { refToString(it,current,master,pulled,true,level) }

        def result = new ArrayList<String>(branches.size() + tags.size())
        if (headRef)
            result.add(refToString(headRef,current,master,pulled,false,level))
        result.addAll(branches)
        result.addAll(tags)
        return result
    }

    @Override
    Map getBranchesAndTags(boolean checkForUpdates) {
        final result = [:]
        if( !hasBareRepo() ){
            return result
        }
        final master = getDefaultBranch()
        final commits = listCommits()
        final branches = []
        final tags = []
        final pulled = new LinkedList<Map>()

        Map<String, Ref> remote = checkForUpdates ? lsRemote() : null
        getBranchList()
                .findAll { it.name.startsWith('refs/heads/') || isRemoteBranch(it) }
                .unique { shortenRefName(it.name) }
                .each {
                    final map = refToMap(it,remote)
                    if (map.commitId as String in commits) pulled << map
                    branches << map
                }

        remote = checkForUpdates ? lsRemote(true) : null
        getTagList()
                .findAll  { it.name.startsWith('refs/tags/') }
                .each {
                    final map = refToMap(it,remote)
                    if (map.commitId as String in commits) pulled << map
                    tags << map
                }

        result.master = master      // master branch name
        result.pulled = pulled.collect { it.name }  // collection of pulled revisions
        result.branches = branches  // collection of branches
        result.tags = tags          // collect of tags
        return result
    }

    @Override
    protected Ref getPeeledRef(Ref ref){
        return bareGitWithLegacyFallback.getRepository().getRefDatabase().peel(ref)
    }

    @Override
    protected List<Ref> getBranchList() {
        bareGitWithLegacyFallback.branchList().setListMode(ListBranchCommand.ListMode.ALL).call()
    }

    @Override
    protected List<Ref> getTagList() {
        bareGitWithLegacyFallback.tagList().call()
    }

    protected String refToString(Ref ref, String current, String master, List<String> pulled, boolean tag, int level ) {

        def result = new StringBuilder()
        def name = shortenRefName(ref.name)
        def peel = getPeeledRef(ref)
        def obj = peel.getPeeledObjectId() ?: peel.getObjectId()
        result << (obj.name in pulled ? 'P' : ' ')

        if( level ) {
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

    @Override
    protected Map<String,Ref> lsRemote(boolean tags = false){
        final cmd = bareGitWithLegacyFallback.lsRemote().setTags(tags)
        if( provider.hasCredentials() )
            cmd.setCredentialsProvider( provider.getGitCredentials() )
        return cmd.callAsMap()
    }

    @Override
    protected String getGitConfigRemoteUrl() {
        if( !bareRepo.exists() ) {
            return null
        }

        final gitConfig = localGitConfig
        if( !gitConfig.exists() ) {
            return null
        }

        final iniFile = new IniFile().load(gitConfig)
        final url = iniFile.getString("remote \"origin\"", "url")
        log.debug "Git config: $gitConfig; url: $url"
        return url
    }

    boolean hasLegacyRepo(){
        return new File(localRootPath, '.git').exists()
    }
}
