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

import static nextflow.util.StringUtils.redactUrlPassword

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.config.Manifest
import nextflow.exception.AbortOperationException
import nextflow.file.FileMutex
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.api.ListBranchCommand
import org.eclipse.jgit.api.errors.RefNotFoundException
import org.eclipse.jgit.errors.RepositoryNotFoundException
import org.eclipse.jgit.internal.storage.file.FileRepository
import org.eclipse.jgit.lib.Constants
import org.eclipse.jgit.lib.ObjectId
import org.eclipse.jgit.lib.Ref
import org.eclipse.jgit.lib.RefUpdate
import org.eclipse.jgit.lib.Repository
import org.eclipse.jgit.lib.SubmoduleConfig
import org.eclipse.jgit.storage.file.FileRepositoryBuilder
import org.eclipse.jgit.transport.FetchResult
import org.eclipse.jgit.transport.RefSpec
import org.eclipse.jgit.transport.TrackingRefUpdate
import org.eclipse.jgit.transport.URIish

/**
 * Local pipeline repository that uses a bare repository with shared clones.
 * This approach allows multiple revisions to coexist efficiently by sharing objects
 * through a bare repository and creating lightweight clones for each commit.
 * The bare repository is stored in `$NXF_ASSETS/.repos/{project}/bare/` and each
 * commit is cloned into `$NXF_ASSETS/.repos/{project}/clones/{commitId}/`.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class MultiRevisionRepository {

    static final String REPOS_SUBDIR = '.repos'
    static final String BARE_REPO = 'bare'
    static final String REVISION_SUBDIR = 'clones'

    private final File root = AssetManager.root

    private String project

    private RepositoryProvider provider

    /**
     * The revision (branch, tag, or commit SHA) to work with
     */
    private String revision

    /**
     * Local directory used in place of the managed clones (e.g. when running a local pipeline directory)
     */
    private File localPathOverride

    private File projectPath
    private File commitPath
    private File bareRepo
    private File revisionSubdir
    private Git _bareGit
    private Git _commitGit
    private Git _localGit

    MultiRevisionRepository(String project, String revision = null) {
        this.project = project
        if( project ) {
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

    void setProvider(RepositoryProvider provider) {
        this.provider = provider
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
            def rev = getBareGit()
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
        final updateRevision = revision ?: getDefaultBranch()
        log.debug "Fetching (updating) bare repo for ${project} [revision: $updateRevision]"
        final fetch = getBareGit().fetch().setRefSpecs(refSpecForName(updateRevision))
        if( provider.hasCredentials() ) {
            fetch.setCredentialsProvider(provider.getGitCredentials())
        }
        verifyFetchResult(fetch.call(), updateRevision, "bare repo")
    }

    private void createBareRepo(Manifest manifest) {
        assert provider
        // This check is required in case of two nextflow instances were doing the same shared clone at the same time.
        // If commitPath exists the previous nextflow instance created the shared clone successfully.
        if( hasBareRepo() )
            return

        try {
            final cloneURL = provider.getCloneUrl()
            log.debug "Pulling bare repo for ${project} -- Using remote clone url: ${redactUrlPassword(cloneURL)}"
            def bare = Git.cloneRepository()
            if( provider.hasCredentials() )
                bare.setCredentialsProvider(provider.getGitCredentials())

            bare
                .setBare(true)
                .setURI(cloneURL)
                .setGitDir(this.bareRepo)
                .setCloneSubmodules(manifest.recurseSubmodules)
                .call()
                .close()
        } catch( Throwable t ) {
            // If there is an error creating the bare repo, remove the bare repo path to avoid incorrect repo clones.
            bareRepo.deleteDir()
            throw t
        }
    }

    /**
     * Download a revision into a shared clone, creating or updating the bare repository as needed.
     *
     * @param revision The revision to download (branch, tag, or commit SHA)
     * @return Status message describing what was done
     */
    String download(String revision, Manifest manifest) {

        // Update revision if specified
        if( revision )
            setRevision(revision)

        // get local copy of bare repository if not exists and fetch revision if specified
        checkBareRepo(manifest)

        // Try to get the commit Id for the revision and abort if not found.
        this.revision ?= getDefaultBranch()
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
        return "downloaded from ${redactUrlPassword(provider.getCloneUrl())}"
    }

    /**
     * Switch to a revision that is already downloaded, by pointing the local path at its clone.
     *
     * @param revision
     * @throws RefNotFoundException
     */
    void tryCheckout(String revision) throws RefNotFoundException {
        if( !hasBareRepo() ) {
            throw new RefNotFoundException("Unknown repository")
        }
        // get the commit ID for the revision. If not specified try to checkout to default branch
        final downloadRevision = revision ?: getDefaultBranch()
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

    Ref checkoutRemoteBranch(String revision, Manifest manifest) {
        download(revision, manifest)
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
            try (FileRepository repo = new FileRepositoryBuilder()
                .setGitDir(new File(commitPath, ".git"))
                .addAlternateObjectDirectory(bareObjectsDir)
                .build() as FileRepository) {
                repo.create()

                // Write alternates file (not done by repo create).
                new File(repo.getObjectsDirectory(), "info/alternates").write(bareObjectsDir.absolutePath)

                // Configure remote pointing to the cache repo
                repo.getConfig().setString("remote", "origin", "url", cloneURL)
                repo.getConfig().save()
            }

            final fetch = getCommitGit()
                .fetch()
                .setRefSpecs(refSpecForName(downloadRevision))
            if( recurseSubmodules )
                fetch.setRecurseSubmodules(SubmoduleConfig.FetchRecurseSubmodulesMode.YES)
            verifyFetchResult(fetch.call(), downloadRevision, "shared clone")
            getCommitGit().checkout().setName(downloadRevision).call()

        } catch( Throwable t ) {
            // If there is an error creating the shared clone, remove the local path to avoid incorrect clones.
            commitPath.deleteDir()
            throw t
        }

        // return status message
        return "downloaded from ${redactUrlPassword(provider.getCloneUrl())}"
    }

    private RefSpec refSpecForName(String revision) {
        // First, check if it's a local branch
        final branchName = "refs/heads/$revision".toString()
        final branch = getBareGit().getRepository().findRef(branchName)
        if( branch != null ) {
            // force-update: branches are mutable and may be force-pushed remotely
            return new RefSpec("+$branchName:$branchName")
        }

        // Check if it's a local tag
        final tagName = "refs/tags/$revision".toString()
        final tag = getBareGit().getRepository().findRef(tagName)
        if( tag != null ) {
            return new RefSpec("$tagName:$tagName")
        }

        // Not found locally - check remote refs
        final remoteRefs = lsRemote(false)

        // Is it a remote branch?
        if( remoteRefs.containsKey(branchName) ) {
            return new RefSpec("+$branchName:$branchName")
        }

        // Is it a remote tag?
        if( remoteRefs.containsKey(tagName) ) {
            return new RefSpec("$tagName:$tagName")
        }

        // Assume it's a commit SHA
        return new RefSpec("$revision:$tagName")
    }

    @PackageScope
    void verifyFetchResult(FetchResult result, String revision, String repoLabel) {
        for( TrackingRefUpdate update : result.getTrackingRefUpdates() ) {
            final r = update.result
            final detail = "${project} [revision: $revision]: ${update.localName} <- ${update.remoteName} => $r (old=${update.oldObjectId?.name()}, new=${update.newObjectId?.name()})"
            if( r == RefUpdate.Result.FORCED ) {
                log.warn "Force-push detected on $repoLabel ($detail); resetting local copy to match remote"
            }
            else if( r == RefUpdate.Result.NEW || r == RefUpdate.Result.FAST_FORWARD || r == RefUpdate.Result.NO_CHANGE ) {
                log.debug "Fetch update on $repoLabel: $detail"
            }
            else {
                throw new AbortOperationException("Unable to update revision '${revision}' in $repoLabel for ${project}: ${update.localName} => $r")
            }
        }
    }

    /**
     * @return The local path of the current revision, or {@code null} if it has not been downloaded
     */
    /**
     * @return The local directory of the current revision. If the revision is not downloaded,
     *         returns a path under the clones directory that does not exist (clones are named by commit id)
     */
    File getLocalPath() {
        if( localPathOverride )
            return localPathOverride
        if( commitPath )
            return commitPath
        return revisionSubdir ? new File(revisionSubdir, revision ?: Constants.HEAD) : null
    }

    protected Git getBareGit() {
        assert bareRepo
        if( !_bareGit ) {
            _bareGit = Git.open(this.bareRepo)
            removeBitbucketToken(_bareGit, provider)
        }
        return _bareGit
    }

    protected Git getCommitGit() {
        if( !_commitGit ) {
            _commitGit = Git.open(commitPath)
        }
        return _commitGit
    }

    private Git getLocalGit() {
        if( !new File(localPathOverride, '.git').exists() )
            return null
        if( !_localGit )
            _localGit = Git.open(localPathOverride)
        return _localGit
    }

    private Git getBareGitOrNull() {
        return hasBareRepo() ? getBareGit() : null
    }

    Git getGit() {
        if( localPathOverride )
            return getLocalGit()
        if( commitPath && commitPath.exists() )
            return getCommitGit()
        return getBareGitOrNull()
    }

    boolean isClean() {
        if( !commitPath?.exists() )
            return true
        try {
            getCommitGit().status().call().isClean()
        }
        catch( RepositoryNotFoundException e ) {
            return true
        }
    }

    /**
     * @return The default branch of the remote repository, as recorded by the HEAD of the bare repository
     */
    String getRemoteDefaultBranch() {
        return getDefaultBranchFromBareRepo()
    }

    String getDefaultBranch() {
        return getRemoteDefaultBranch() ?: Const.DEFAULT_BRANCH
    }

    List<Ref> getBranchList() {
        getBareGitOrNull()?.branchList()?.setListMode(ListBranchCommand.ListMode.ALL)?.call() ?: []
    }

    List<Ref> getTagList() {
        getBareGitOrNull()?.tagList()?.call() ?: []
    }

    Map<String, Ref> lsRemote(boolean tags) {
        final cmd = getBareGitOrNull()?.lsRemote()?.setTags(tags)
        if( provider?.hasCredentials() )
            cmd?.setCredentialsProvider(provider.getGitCredentials())
        return cmd?.callAsMap() ?: [:]
    }

    Ref peel(Ref ref) {
        return getBareGitOrNull()?.getRepository()?.getRefDatabase()?.peel(ref)
    }

    File getLocalGitConfig() {
        if( localPathOverride )
            return new File(localPathOverride, '.git/config')
        return hasBareRepo() ? new File(this.bareRepo, 'config') : null
    }

    List<String> listDownloadedCommits() {
        def result = new LinkedList()
        if( !getRevisionSubdir()?.exists() )
            return result

        getRevisionSubdir().eachDir { File it -> result << it.getName().toString() }

        return result
    }

    /**
     * Drop local copy of a repository. If revision is specified, only removes the specified revision
     */
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
                throw new AbortOperationException("Unable to delete project `${project}` -- Check access permissions for path: ${commitPath}")
            return
        }

        throw new AbortOperationException("Revision $revision for project `${project}` contains uncommitted changes -- won't drop it")
    }

    void close() {
        if( _bareGit ) {
            _bareGit.close()
            _bareGit = null
        }
        if( _commitGit ) {
            _commitGit.close()
            _commitGit = null
        }
        if( _localGit ) {
            _localGit.close()
            _localGit = null
        }
    }

    void setLocalPath(File file) {
        localPathOverride = file
    }

    void setProject(String project) {
        this.project = project
        this.projectPath = new File(root, REPOS_SUBDIR + '/' + project)
        this.bareRepo = new File(projectPath, BARE_REPO)
        this.revisionSubdir = new File(projectPath, REVISION_SUBDIR)
    }

    File getProjectPath() {
        return projectPath
    }

    String getCurrentRevision() {
        Ref head = findHeadRef()
        if( !head )
            return '(unknown)'

        if( head.isSymbolic() )
            return Repository.shortenRefName(head.getTarget().getName())

        if( !head.getObjectId() )
            return '(unknown)'

        // try to resolve the current object id to a tag name
        def name = resolveTagNameByObjectId(head.objectId)
        return name ? Repository.shortenRefName(name) : head.objectId.name()
    }

    AssetManager.RevisionInfo getCurrentRevisionAndName() {
        Ref head = findHeadRef()
        if( !head )
            return null

        if( head.isSymbolic() ) {
            return new AssetManager.RevisionInfo(
                head.objectId.name(),
                Repository.shortenRefName(head.getTarget().getName()),
                AssetManager.RevisionInfo.Type.BRANCH
            )
        }

        if( !head.getObjectId() )
            return null

        final name = resolveTagNameByObjectId(head.objectId)
        if( name ) {
            return new AssetManager.RevisionInfo(head.objectId.name(), Repository.shortenRefName(name), AssetManager.RevisionInfo.Type.TAG)
        } else {
            return new AssetManager.RevisionInfo(head.objectId.name(), null, AssetManager.RevisionInfo.Type.COMMIT)
        }
    }

    protected Ref findHeadRef() {
        getGit()?.getRepository()?.findRef(Constants.HEAD)
    }

    protected String resolveTagNameByObjectId(ObjectId objectId) {
        Collection<Ref> tags = getGit()?.getRepository()?.getRefDatabase()?.getRefsByPrefix(Constants.R_TAGS)
        return tags?.find { it.objectId == objectId || it.peeledObjectId == objectId }?.name
    }

    /**
     * Remove the Bitbucket API token from the remote URL of a local repository, if present.
     *
     * Previous versions stored the Bitbucket API token in the clone URL, see
     * https://github.com/nextflow-io/nextflow/issues/7868
     */
    @PackageScope
    static void removeBitbucketToken(Git git, RepositoryProvider provider) {
        if( provider !instanceof BitbucketRepositoryProvider || !provider.getToken() )
            return
        try {
            final config = git.getRepository().getConfig()
            final url = config.getString('remote', 'origin', 'url')
            if( !url )
                return
            final uri = new URIish(url)
            if( !uri.getPass() )
                return
            config.setString('remote', 'origin', 'url', uri.setPass(null).toString())
            config.save()
        }
        catch( Exception e ) {
            log.debug "Unable to remove Bitbucket token from remote url -- cause: ${e.message}"
        }
    }
}
