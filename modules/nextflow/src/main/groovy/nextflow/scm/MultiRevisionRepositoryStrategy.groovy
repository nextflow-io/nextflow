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
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.config.Manifest
import nextflow.exception.AbortOperationException
import nextflow.file.FileMutex
import nextflow.util.IniFile
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

    static final String BARE_REPO = '.nextflow/bare_repo'
    static final String REVISION_SUBDIR = '.nextflow/commits'

    /**
     * The revision (branch, tag, or commit SHA) to work with
     */
    private String revision

    /**
     * Path to the commit-specific clone
     */
    private File commitPath
    private File bareRepo
    private File revisionSubdir
    private Git _bareGit
    private Git _commitGit

    MultiRevisionRepositoryStrategy(RepositoryContext context, String revision = null) {
        super(context)
        this.bareRepo = new File(context.root, context.project + '/' + BARE_REPO)
        this.revisionSubdir = new File(context.root, context.project + '/' + REVISION_SUBDIR)
        if (revision)
            setRevision(revision)
    }

    @PackageScope
    File getBareRepo() { this.bareRepo }

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
        if( !commitId )
            return
        final oldCommitPath = this.commitPath
        if( oldCommitPath && _commitGit ) {
            _commitGit.close()
            _commitGit = null
        }
        this.commitPath = new File(getRevisionSubdir(), commitId)
    }

    @PackageScope
    String revisionToCommitWithBareRepo(String revision) {
        String commitId = null

        if( getBareRepo().exists() ) {
            def rev = Git.open(getBareRepo())
                .getRepository()
                .resolve(revision ?: Constants.HEAD)
            if( rev )
                commitId = rev.getName()
        }

        return commitId
    }

    protected String getDefaultBranchFromBareRepo() {
        if( !getBareRepo().exists() ) {
            log.debug "Bare repo ${getBareRepo()} doesn't exist"
            return null
        }
        def head = getDefaultBranchRef()
        if( head ) {
            Repository.shortenRefName(head.getName())
        } else {
            log.debug "Unable to determine default branch from bare repo ${getBareRepo()}"
            return null
        }
    }

    protected Ref getDefaultBranchRef() {
        return getBareGit().getRepository().findRef(Constants.HEAD).getTarget()
    }

    protected void checkBareRepo(Manifest manifest) {
        /*
         * if the bare repository of the pipeline does not exists locally pull it from the remote repo
         */
        if( !getBareRepo().exists() ) {
            getBareRepo().parentFile.mkdirs()
            // Use a file mutex to prevent concurrent clones of the same commit.
            final file = new File(getBareRepo().parentFile, ".${getBareRepo().name}.lock")
            final wait = "Another Nextflow instance is creating the bare repo for ${context.project} -- please wait till it completes"
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
        log.debug "Fetching (updating) bare repo for ${context.project} [revision: $updateRevision]"
        getBareGit().fetch().setRefSpecs(refSpecForName(updateRevision)).call()
    }

    private void createBareRepo(Manifest manifest) {
        // This check is required in case of two nextflow instances were doing the same shared clone at the same time.
        // If commitPath exists the previous nextflow instance created the shared clone successfully.
        if( getBareRepo().exists() )
            return

        try {
            final cloneURL = getGitRepositoryUrl()
            log.debug "Pulling bare repo for ${context.project} -- Using remote clone url: ${cloneURL}"
            def bare = Git.cloneRepository()
            if( context.provider.hasCredentials() )
                bare.setCredentialsProvider(context.provider.getGitCredentials())

            bare
                .setBare(true)
                .setURI(cloneURL)
                .setGitDir(getBareRepo())
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
            throw new AbortOperationException("No commit found for revision ${this.revision} in project ${context.project}")
        }
        updateCommitDir(commitId)

        if( commitPath.exists() )
            return "already downloaded from ${getGitRepositoryUrl()}"

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
        if( !bareRepo.exists() ) {
            throw new RefNotFoundException("Unknown repository")
        }
        // get the commit ID for the revision. If not specified try to checkout to default branch
        final downloadRevision = revision ?: getDefaultBranch(manifest)
        final commitId = revisionToCommitWithBareRepo(downloadRevision)
        if( !commitId ) {
            throw new RefNotFoundException("No commit found for revision ${downloadRevision} in project ${context.project}")
        }
        if( new File(revisionSubdir, commitId).exists() ) {
            setRevision(downloadRevision)
        } else {
            throw new RefNotFoundException("Revision ${downloadRevision} not locally checked out.")
        }
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
            return "already downloaded"

        commitPath.parentFile.mkdirs()
        try {
            final cloneURL = getBareRepo().toString()
            log.debug "Pulling ${context.project} -- Using remote clone url: ${cloneURL}"

            // clone it, but don't specify a revision - jgit will checkout the default branch
            File bareObjectsDir = new File(getBareRepo(), "objects")
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
        // Is it a local branch?
        Ref branch = getBareGit().getRepository().findRef("refs/heads/" + revision)
        if( branch != null ) {
            return new RefSpec("refs/heads/" + revision + ":refs/heads/" + revision)
        }

        // Is it a tag?
        Ref tag = getBareGit().getRepository().findRef("refs/tags/" + revision)
        if( tag != null ) {
            return new RefSpec("refs/tags/" + revision + ":refs/tags/" + revision)
        }

        // It is a commit
        return new RefSpec(revision + ":refs/tags/" + revision)
    }

    @Override
    File getLocalPath() {
        return this.commitPath ?: context.localRootPath
    }

    protected Git getBareGit() {
        if( !_bareGit ) {
            _bareGit = Git.open(getBareRepo())
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
        if( getBareRepo().exists() )
            return getBareGit()
        // Fallback to legacy if it exists
        if( context.localRootPath.exists() && new File(context.localRootPath, '.git').exists() )
            return Git.open(context.localRootPath)
        return null
    }

    private Git getBareGitWithLegacyFallback() {
        if( getBareRepo().exists() )
            return getBareGit()
        // Fallback to legacy
        if( context.localRootPath.exists() && new File(context.localRootPath, '.git').exists() )
            return Git.open(context.localRootPath)
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
        if( context.provider.hasCredentials() )
            cmd?.setCredentialsProvider(context.provider.getGitCredentials())
        return cmd?.callAsMap() ?: [:]
    }

    @Override
    Ref peel(Ref ref) {
        return getBareGitWithLegacyFallback()?.getRepository()?.getRefDatabase()?.peel(ref)
    }

    @Override
    File getLocalGitConfig() {
        getBareRepo().exists() ? new File(getBareRepo(), 'config')
            : (context.localRootPath ? new File(context.localRootPath, '.git/config') : null)
    }

    @Override
    List<String> listDownloadedCommits() {
        def result = new LinkedList()
        if( !context.localRootPath.exists() )
            return result
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
            context.localRootPath.deleteDir()
        }
    }

    private void dropRevision(String revision, boolean force) {
        assert revision
        setRevision(revision)
        if( force || isClean() ) {
            close()
            if( !localPath.deleteDir() )
                throw new AbortOperationException("Unable to delete project `${context.project}` -- Check access permissions for path: ${localPath}")
            return
        }

        throw new AbortOperationException("Revision $revision for project `${context.project}` contains uncommitted changes -- won't drop it")
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
}
