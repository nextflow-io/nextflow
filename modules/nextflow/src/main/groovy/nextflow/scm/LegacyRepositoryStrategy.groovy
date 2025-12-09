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
import groovy.util.logging.Slf4j
import nextflow.config.Manifest
import nextflow.exception.AbortOperationException
import org.eclipse.jgit.api.CreateBranchCommand
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.api.ListBranchCommand
import org.eclipse.jgit.api.MergeResult
import org.eclipse.jgit.api.errors.RefNotFoundException
import org.eclipse.jgit.errors.RepositoryNotFoundException
import org.eclipse.jgit.lib.Constants
import org.eclipse.jgit.lib.Ref
import org.eclipse.jgit.lib.Repository
import org.eclipse.jgit.lib.SubmoduleConfig.FetchRecurseSubmodulesMode

/**
 * Legacy repository strategy that directly clones repositories to the local path.
 * This is the traditional approach where each project gets a full git clone.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class LegacyRepositoryStrategy extends AbstractRepositoryStrategy {

    private Git _git

    LegacyRepositoryStrategy(RepositoryContext context) {
        super(context)
    }

    @Override
    String download(String revision, Integer deep, Manifest manifest) {
        /*
         * if the pipeline already exists locally pull it from the remote repo
         */
        if( !getLocalPath().exists() ) {
            getLocalPath().parentFile.mkdirs()

            final cloneURL = getGitRepositoryUrl()
            log.debug "Pulling ${context.project} -- Using remote clone url: ${cloneURL}"

            // clone it, but don't specify a revision - jgit will checkout the default branch
            def clone = Git.cloneRepository()
            if( context.provider.hasCredentials() )
                clone.setCredentialsProvider( context.provider.getGitCredentials() )

            clone
                .setURI(cloneURL)
                .setDirectory(getLocalPath())
                .setCloneSubmodules(manifest.recurseSubmodules)
            if( deep )
                clone.setDepth(deep)
            clone.call()

            // git cli would automatically create a 'refs/remotes/origin/HEAD' symbolic ref pointing at the remote's
            // default branch. jgit doesn't do this, but since it automatically checked out the default branch on clone
            // we can create the symbolic ref ourselves using the current head
            def head = getGit().getRepository().findRef(Constants.HEAD)
            if( head ) {
                def headName = head.isSymbolic()
                    ? Repository.shortenRefName(head.getTarget().getName())
                    : head.getName()

                getGit().repository.getRefDatabase()
                    .newUpdate(REMOTE_DEFAULT_HEAD, true)
                    .link(REMOTE_REFS_ROOT + headName)
            } else {
                log.debug "Unable to determine default branch of repo ${cloneURL}, symbolic ref not created"
            }

            // now the default branch is recorded in the repo, explicitly checkout the revision (if specified).
            // this also allows 'revision' to be a SHA commit id, which isn't supported by the clone command
            if( revision ) {
                try { getGit().checkout() .setName(revision) .call() }
                catch ( RefNotFoundException e ) { checkoutRemoteBranch(revision, manifest) }
            }

            // return status message
            return "downloaded from ${cloneURL}"
        }

        log.debug "Pull pipeline ${context.project}  -- Using local path: ${getLocalPath()}"

        // verify that is clean
        if( !isClean() )
            throw new AbortOperationException("${context.project} contains uncommitted changes -- cannot pull from repository")

        if( revision && revision != getCurrentRevision() ) {
            /*
             * check out a revision before the pull operation
             */
            try {
                getGit().checkout() .setName(revision) .call()
            }
            /*
             * If the specified revision does not exist
             * Try to checkout it from a remote branch and return
             */
            catch ( RefNotFoundException e ) {
                final ref = checkoutRemoteBranch(revision, manifest)
                final commitId = ref?.getObjectId()
                return commitId
                    ? "checked out at ${commitId.name()}"
                    : "checked out revision ${revision}"
            }
        }

        def pull = getGit().pull()
        def revInfo = getCurrentRevisionAndName()

        if ( revInfo.type == AssetManager.RevisionInfo.Type.COMMIT ) {
            log.debug("Repo appears to be checked out to a commit hash, but not a TAG, so we will assume the repo is already up to date and NOT pull it!")
            return MergeResult.MergeStatus.ALREADY_UP_TO_DATE.toString()
        }

        if ( revInfo.type == AssetManager.RevisionInfo.Type.TAG ) {
            pull.setRemoteBranchName( "refs/tags/" + revInfo.name )
        }

        if( context.provider.hasCredentials() )
            pull.setCredentialsProvider( context.provider.getGitCredentials() )

        if( manifest.recurseSubmodules ) {
            pull.setRecurseSubmodules(FetchRecurseSubmodulesMode.YES)
        }
        def result = pull.call()
        if(!result.isSuccessful())
            throw new AbortOperationException("Cannot pull project `${context.project}` -- ${result.toString()}")

        return result?.mergeResult?.mergeStatus?.toString()
    }

    @Override
    void tryCheckout(String revision, Manifest manifest) throws RefNotFoundException {
        assert localPath

        def current = getCurrentRevision()
        if( current != getDefaultBranch(manifest) ) {
            if( !revision ) {
                throw new AbortOperationException("Project `$context.project` is currently stuck on revision: $current -- you need to explicitly specify a revision with the option `-r` in order to use it")
            }
        }
        if( !revision || revision == current ) {
            // nothing to do
            return
        }

        // verify that is clean
        if( !isClean() )
            throw new AbortOperationException("Project `$context.project` contains uncommitted changes -- Cannot switch to revision: $revision")


        git.checkout().setName(revision) .call()
    }

    @Override
    Ref checkoutRemoteBranch(String revision, Manifest manifest) {
        try {
            def fetch = git.fetch()
            if(context.provider.hasCredentials()) {
                fetch.setCredentialsProvider( context.provider.getGitCredentials() )
            }
            if( manifest.recurseSubmodules ) {
                fetch.setRecurseSubmodules(FetchRecurseSubmodulesMode.YES)
            }
            fetch.call()

            try {
                return git.checkout()
                    .setCreateBranch(true)
                    .setName(revision)
                    .setUpstreamMode(CreateBranchCommand.SetupUpstreamMode.TRACK)
                    .setStartPoint("origin/" + revision)
                    .call()
            }
            catch (RefNotFoundException e) {
                return git.checkout() .setName(revision) .call()
            }
        }
        catch (RefNotFoundException e) {
            throw new AbortOperationException("Cannot find revision `$revision` -- Make sure that it exists in the remote repository `$gitRepositoryUrl`", e)
        }
    }

    @Override
    File getLocalPath() {
        return context.localRootPath
    }

    @Override
    Git getGit() {
        if( !_git ) {
            _git = Git.open(getLocalPath())
        }
        return _git
    }

    @Override
    boolean isClean() {
        try {
            getGit().status().call().isClean()
        }
        catch( RepositoryNotFoundException e ) {
            return true
        }
    }

    @Override
    String getRemoteDefaultBranch() {
        Ref remoteHead = git.getRepository().findRef(REMOTE_DEFAULT_HEAD)
        return remoteHead?.getTarget()?.getName()?.substring(REMOTE_REFS_ROOT.length())
    }

    @Override
    List<Ref> getBranchList() {
        getGit().branchList().setListMode(ListBranchCommand.ListMode.ALL).call()
    }

    @Override
    List<Ref> getTagList() {
        getGit().tagList().call()
    }

    @Override
    Map<String, Ref> lsRemote(boolean tags) {
        final cmd = getGit().lsRemote().setTags(tags)
        if( context.provider.hasCredentials() )
            cmd.setCredentialsProvider( context.provider.getGitCredentials() )
        return cmd.callAsMap()
    }

    @Override
    Ref peel(Ref ref) {
        return getGit().getRepository().getRefDatabase().peel(ref)
    }

    @Override
    File getLocalGitConfig() {
        getLocalPath() ? new File(getLocalPath(),'.git/config') : null
    }

    @Override
    List<String> listDownloadedCommits() {
        if( !context.localRootPath.exists() )
            return []
        return [findHeadRef().objectId.getName()]
    }

    @Override
    void drop(String revision, boolean force) {
        if (!localPath.exists())
            throw new AbortOperationException("No match found for: ${context.project}")

        if (revision)
            throw new AbortOperationException("Not able to remove a revision for a Legacy Repo. Use all option to remove the local repository.")

        if (force || isClean()) {
            close()
            if( !localPath.deleteDir() )
                throw new AbortOperationException("Unable to delete project `${context.project}` -- Check access permissions for path: ${localPath}")
            return
        }
        throw new AbortOperationException("Local project repository contains uncommitted changes -- won't drop it")
    }

    @Override
    void close() {
        if( _git ) {
            _git.close()
            _git = null
        }
    }
}
