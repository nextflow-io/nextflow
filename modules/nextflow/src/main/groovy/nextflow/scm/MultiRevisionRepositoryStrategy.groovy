/*
 * Copyright 2013-2025, Seqera Labs
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

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.config.Manifest
import nextflow.exception.AbortOperationException
import nextflow.file.FileMutex
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.api.ListBranchCommand
import org.eclipse.jgit.api.errors.RefNotFoundException
import org.eclipse.jgit.errors.RepositoryNotFoundException
import org.eclipse.jgit.internal.storage.file.FileRepository
import org.eclipse.jgit.lib.Constants
import org.eclipse.jgit.lib.Ref
import org.eclipse.jgit.lib.Repository
import org.eclipse.jgit.lib.SubmoduleConfig
import org.eclipse.jgit.storage.file.FileRepositoryBuilder
import org.eclipse.jgit.transport.RefSpec

/**
 * Multi-revision repository strategy that uses a bare repository with shared clones.
 * This approach allows multiple revisions to coexist efficiently by sharing objects
 * through a bare repository and creating lightweight clones for each commit.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class MultiRevisionRepositoryStrategy extends AbstractRepositoryStrategy {

    static final String REPOS_SUBDIR = '.repos'
    static final String BARE_REPO = 'bare'
    static final String REVISION_SUBDIR = 'clones'

    /**
     * The revision (branch, tag, or commit SHA) to work with
     */
    private String revision

    /**
     * Path to the commit-specific clone
     */
    private File legacyRepoPath
    private File projectPath
    private File commitPath
    private File bareRepo
    private File revisionSubdir
    private Git _bareGit
    private Git _commitGit

    MultiRevisionRepositoryStrategy(String project, String revision = null) {
        super(project)
        if( project ) {
            this.legacyRepoPath = new File(root, project)
            this.projectPath = new File(root, REPOS_SUBDIR + '/' + project)
            this.bareRepo = new File(projectPath, BARE_REPO)
            this.revisionSubdir = new File(projectPath, REVISION_SUBDIR)
        }
        if( revision )
            setRevision(revision)
        else if( hasBareRepo() ) {
            setRevision(getDefaultBranchFromBareRepo())
        }

    }
    /**
     * Check if a project exists and has the bare repo for multi-revision
     * @param root
     * @param project
     * @return
     */
    static boolean checkProject(File root, String project) {
        return new File(root, REPOS_SUBDIR + '/' + project + '/' + BARE_REPO).exists()
    }

    @PackageScope
    File getBareRepo() { this.bareRepo }

    boolean hasBareRepo() {
        return this.bareRepo && this.bareRepo.exists()
    }

    @PackageScope
    File getRevisionSubdir() { this.revisionSubdir }

    @PackageScope
    void setRevision(String revision) {
        assert revision
        this.revision = revision
        updateCommitDir(revisionToCommitWithBareRepo(revision))
    }

    String getRevision() {
        return revision
    }

    @PackageScope
    void updateCommitDir(String commitId) {
        final oldCommitPath = this.commitPath
        if( oldCommitPath && _commitGit ) {
            _commitGit.close()
            _commitGit = null
        }
        this.commitPath = commitId ? new File(getRevisionSubdir(), commitId) : null
    }

    @PackageScope
    String revisionToCommitWithBareRepo(String revision) {
        String commitId = null

        if( hasBareRepo() ) {
            def rev = Git.open(this.bareRepo)
                .getRepository()
                .resolve(revision ?: Constants.HEAD)
            if( rev )
                commitId = rev.getName()
        }

        return commitId
    }

    protected String getDefaultBranchFromBareRepo() {
        if( !hasBareRepo() ) {
            log.debug "Bare repo ${this.bareRepo} doesn't exist"
            return null
        }
        def head = getDefaultBranchRef()
        if( head ) {
            Repository.shortenRefName(head.getName())
        } else {
            log.debug "Unable to determine default branch from bare repo ${this.bareRepo}"
            return null
        }
    }

    protected Ref getDefaultBranchRef() {
        return getBareGit().getRepository().findRef(Constants.HEAD).getTarget()
    }

    protected void checkBareRepo(Manifest manifest) {
        assert bareRepo
        /*
         * if the bare repository of the pipeline does not exists locally pull it from the remote repo
         */
        if( !hasBareRepo() ) {
            if( !bareRepo.parentFile.exists() )
                this.bareRepo.parentFile.mkdirs()
            // Use a file mutex to prevent concurrent clones of the same commit.
            final file = new File(this.bareRepo.parentFile, ".${this.bareRepo.name}.lock")
            final wait = "Another Nextflow instance is creating the bare repo for ${project} -- please wait till it completes"
            final err = "Unable to acquire exclusive lock after 60s on file: $file"

            final mutex = new FileMutex(target: file, timeout: '60s', waitMessage: wait, errorMessage: err)
            try {
                mutex.lock { createBareRepo(manifest) }
            }
            finally {
                file.delete()
            }
        }
        final updateRevision = revision ?: getDefaultBranch(manifest)
        log.debug "Fetching (updating) bare repo for ${project} [revision: $updateRevision]"
        final fetch = getBareGit().fetch().setRefSpecs(refSpecForName(updateRevision))
        if( provider.hasCredentials() ) {
            fetch.setCredentialsProvider(provider.getGitCredentials())
        }
        fetch.call()
    }

    private void createBareRepo(Manifest manifest) {
        assert provider
        // This check is required in case of two nextflow instances were doing the same shared clone at the same time.
        // If commitPath exists the previous nextflow instance created the shared clone successfully.
        if( hasBareRepo() )
            return

        try {
            final cloneURL = provider.getCloneUrl()
            log.debug "Pulling bare repo for ${project} -- Using remote clone url: ${cloneURL}"
            def bare = Git.cloneRepository()
            if( provider.hasCredentials() )
                bare.setCredentialsProvider(provider.getGitCredentials())

            bare
                .setBare(true)
                .setURI(cloneURL)
                .setGitDir(this.bareRepo)
                .setCloneSubmodules(manifest.recurseSubmodules)
                .call()
        } catch( Throwable t ) {
            // If there is an error creating the bare repo, remove the bare repo path to avoid incorrect repo clones.
            bareRepo.deleteDir()
            throw t
        }
    }

    @Override
    String download(String revision, Integer deep, Manifest manifest) {

        // Update revision if specified
        if( revision )
            setRevision(revision)

        // get local copy of bare repository if not exists and fetch revision if specified
        checkBareRepo(manifest)

        // Try to get the commit Id for the revision and abort if not found.
        this.revision ?= getDefaultBranch(manifest)
        final commitId = revisionToCommitWithBareRepo(this.revision)
        if( !commitId ) {
            throw new AbortOperationException("No commit found for revision ${this.revision} in project ${project}")
        }
        updateCommitDir(commitId)

        if( commitPath.exists() )
            return "Already-up-to-date"

        /*
         * if revision does not exists locally pull it from the remote repo
         */
        if( !commitPath.parentFile.exists() )
            commitPath.parentFile.mkdirs()

        // Use a file mutex to prevent concurrent clones of the same commit.
        final file = new File(commitPath.parentFile, ".${commitPath.name}.lock")
        final wait = "Another Nextflow instance is creating clone for $this.revision -- please wait till it completes"
        final err = "Unable to acquire exclusive lock after 60s on file: $file"

        final mutex = new FileMutex(target: file, timeout: '60s', waitMessage: wait, errorMessage: err)
        try {
            mutex.lock { createSharedClone(this.revision, manifest.recurseSubmodules) }
        }
        finally {
            file.delete()
        }
        return "downloaded from ${getGitRepositoryUrl()}"
    }

    /** Checkout in multi-revision is just updating the commit dir. This methods check if bare repo, commit id and commit path exists.
     * If exists updates the revision that updates the localPath reference. If some
     *
     * @param revision
     * @param recurseSubmodules
     * @throws RefNotFoundException
     */
    @Override
    void tryCheckout(String revision, Manifest manifest) throws RefNotFoundException {
        if( !hasBareRepo() ) {
            throw new RefNotFoundException("Unknown repository")
        }
        // get the commit ID for the revision. If not specified try to checkout to default branch
        final downloadRevision = revision ?: getDefaultBranch(manifest)
        final commitId = revisionToCommitWithBareRepo(downloadRevision)
        if( !commitId ) {
            throw new RefNotFoundException("No commit found for revision ${downloadRevision} in project ${project}")
        }
        if( isRevisionLocal(commitId) ) {
            setRevision(downloadRevision)
        } else {
            throw new RefNotFoundException("Revision ${downloadRevision} not locally checked out.")
        }
    }

    private boolean isRevisionLocal(String commitId) {
        this.revisionSubdir && new File(revisionSubdir, commitId).exists()
    }

    @Override
    Ref checkoutRemoteBranch(String revision, Manifest manifest) {
        download(revision, 1, manifest)
        return findHeadRef()
    }

    private String createSharedClone(String downloadRevision, boolean recurseSubmodules) {
        // This check is required in case of two nextflow instances were doing the same shared clone at the same time.
        // If commitPath exists the previous nextflow instance created the shared clone successfully.
        if( commitPath.exists() )
            return "Already-up-to-date"

        commitPath.parentFile.mkdirs()
        try {
            final cloneURL = this.bareRepo.toString()
            log.debug "Pulling ${project} -- Using remote clone url: ${cloneURL}"

            // clone it, but don't specify a revision - jgit will checkout the default branch
            File bareObjectsDir = new File(this.bareRepo, "objects")
            FileRepository repo = new FileRepositoryBuilder()
                .setGitDir(new File(getLocalPath(), ".git"))
                .addAlternateObjectDirectory(bareObjectsDir)
                .build() as FileRepository
            repo.create()

            // Write alternates file (not done by repo create).
            new File(repo.getObjectsDirectory(), "info/alternates").write(bareObjectsDir.absolutePath)

            // Configure remote pointing to the cache repo
            repo.getConfig().setString("remote", "origin", "url", cloneURL)
            repo.getConfig().save()

            final fetch = getCommitGit()
                .fetch()
                .setRefSpecs(refSpecForName(downloadRevision))
            if( recurseSubmodules )
                fetch.setRecurseSubmodules(SubmoduleConfig.FetchRecurseSubmodulesMode.YES)
            fetch.call()
            getCommitGit().checkout().setName(downloadRevision).call()

        } catch( Throwable t ) {
            // If there is an error creating the shared clone, remove the local path to avoid incorrect clones.
            commitPath.deleteDir()
            throw t
        }

        // return status message
        return "downloaded from ${getGitRepositoryUrl()}"
    }

    private RefSpec refSpecForName(String revision) {
        // First, check if it's a local branch
        Ref branch = getBareGit().getRepository().findRef("refs/heads/" + revision)
        if( branch != null ) {
            return new RefSpec("refs/heads/" + revision + ":refs/heads/" + revision)
        }

        // Check if it's a local tag
        Ref tag = getBareGit().getRepository().findRef("refs/tags/" + revision)
        if( tag != null ) {
            return new RefSpec("refs/tags/" + revision + ":refs/tags/" + revision)
        }

        // Not found locally - check remote refs
        final remoteRefs = lsRemote(false)

        // Is it a remote branch?
        if( remoteRefs.containsKey("refs/heads/" + revision) ) {
            return new RefSpec("refs/heads/" + revision + ":refs/heads/" + revision)
        }

        // Is it a remote tag?
        if( remoteRefs.containsKey("refs/tags/" + revision) ) {
            return new RefSpec("refs/tags/" + revision + ":refs/tags/" + revision)
        }

        // Assume it's a commit SHA
        return new RefSpec(revision + ":refs/tags/" + revision)
    }

    @Override
    File getLocalPath() {
        return this.commitPath ?: legacyRepoPath
    }

    protected Git getBareGit() {
        assert bareRepo
        if( !_bareGit ) {
            _bareGit = Git.open(this.bareRepo)
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
    Git getGit() {
        if( commitPath && commitPath.exists() )
            return getCommitGit()
        if( hasBareRepo() )
            return getBareGit()
        if( this.legacyRepoPath && new File(this.legacyRepoPath, '.git').exists() )
            return Git.open(this.legacyRepoPath)
        return null
    }

    private Git getBareGitWithLegacyFallback() {
        if( hasBareRepo() )
            return getBareGit()
        // Fallback to legacy
        if( this.legacyRepoPath && new File(this.legacyRepoPath, '.git').exists() )
            return Git.open(this.legacyRepoPath)
        return null
    }

    @Override
    boolean isClean() {
        try {
            getCommitGit().status().call().isClean()
        }
        catch( RepositoryNotFoundException e ) {
            return true
        }
    }

    @Override
    String getRemoteDefaultBranch() {
        return getDefaultBranchFromBareRepo() ?: findRemoteDefaultBranch()
    }

    String findRemoteDefaultBranch() {
        Ref remoteHead = git.getRepository().findRef(REMOTE_DEFAULT_HEAD)
        return remoteHead?.getTarget()?.getName()?.substring(REMOTE_REFS_ROOT.length())
    }

    @Override
    List<Ref> getBranchList() {
        getBareGitWithLegacyFallback()?.branchList()?.setListMode(ListBranchCommand.ListMode.ALL)?.call() ?: []
    }

    @Override
    List<Ref> getTagList() {
        getBareGitWithLegacyFallback()?.tagList()?.call() ?: []
    }

    @Override
    Map<String, Ref> lsRemote(boolean tags) {
        final cmd = getBareGitWithLegacyFallback()?.lsRemote()?.setTags(tags)
        if( provider?.hasCredentials() )
            cmd?.setCredentialsProvider(provider.getGitCredentials())
        return cmd?.callAsMap() ?: [:]
    }

    @Override
    Ref peel(Ref ref) {
        return getBareGitWithLegacyFallback()?.getRepository()?.getRefDatabase()?.peel(ref)
    }

    @Override
    File getLocalGitConfig() {
        return hasBareRepo() ? new File(this.bareRepo, 'config') : legacyRepoPath ? new File(legacyRepoPath,'.git/config') : null
    }

    @Override
    String getGitRepositoryUrl() {
        return provider.getCloneUrl()
    }

    @Override
    List<String> listDownloadedCommits() {
        def result = new LinkedList()
        if( !getRevisionSubdir().exists() )
            return result

        getRevisionSubdir().eachDir { File it -> result << it.getName().toString() }

        return result
    }

    @Override
    void drop(String revision, boolean force) {
        if( revision ) {
            dropRevision(revision, force)
        } else {
            listDownloadedCommits().each { dropRevision(it, force) }
            bareRepo.parentFile.deleteDir()
        }
    }

    private void dropRevision(String revision, boolean force) {
        assert revision
        setRevision(revision)

        if( !commitPath || !commitPath.exists() ) {
            log.info "No local folder found for revision '$revision' in '$project' -- Nothing to do"
            return
        }

        if( force || isClean() ) {
            close()
            if( !commitPath.deleteDir() )
                throw new AbortOperationException("Unable to delete project `${project}` -- Check access permissions for path: ${localPath}")
            return
        }

        throw new AbortOperationException("Revision $revision for project `${project}` contains uncommitted changes -- won't drop it")
    }

    @Override
    void close() {
        if( _bareGit ) {
            _bareGit.close()
            _bareGit = null
        }
        if( _commitGit ) {
            _commitGit.close()
            _commitGit = null
        }
    }

    @Override
    void setLocalPath(File file) {
        legacyRepoPath = file
    }

    @Override
    void setProject(String project) {
        super.setProject(project)
        this.projectPath = new File(root, REPOS_SUBDIR + '/' + project)
        this.bareRepo = new File(projectPath, BARE_REPO)
        this.revisionSubdir = new File(projectPath, REVISION_SUBDIR)
    }

    @Override
    File getProjectPath() {
        return projectPath
    }
}
