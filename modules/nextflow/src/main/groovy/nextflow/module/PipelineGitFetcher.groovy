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

package nextflow.module

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.api.errors.GitAPIException
import org.eclipse.jgit.api.errors.JGitInternalException
import org.eclipse.jgit.transport.RefSpec
import org.eclipse.jgit.transport.URIish

/**
 * Fetches the exact tree of a pinned commit from a public Git repository.
 *
 * Only public repositories are supported (no credentials). Switch to
 * AssetManager if private registries arrive.
 */
@Slf4j
@CompileStatic
class PipelineGitFetcher {

    /**
     * Fetch the pinned commit, or the revision when the release has no commit id
     *
     * @param keepGit Keep the `.git` folder, with `origin` set to the repository URL
     * @param depth The history depth to fetch, 1 by default
     * @return The commit that was fetched
     */
    String fetch(String repositoryUrl, String revision, String commitId, Path target, boolean keepGit = false, Integer depth = null) {
        if( !revision && !commitId )
            throw new AbortOperationException("Cannot fetch $repositoryUrl: the release has neither a git revision nor a commit id")
        Git git = null
        String head
        try {
            // the revision can be a branch that moved after the release, so the commit takes precedence
            git = commitId
                ? fetchCommit(repositoryUrl, commitId, target, depth ?: 1)
                : Git.cloneRepository().setURI(repositoryUrl).setDirectory(target.toFile()).setBranch(revision).setDepth(depth ?: 1).call()
            head = git.repository.resolve('HEAD')?.name()
        }
        catch( GitAPIException | JGitInternalException e ) {
            final ref = commitId ? "commit '${commitId}'" : "revision '${revision}'"
            throw new AbortOperationException("Cannot fetch ${ref} of $repositoryUrl: ${e.message}", e)
        }
        finally {
            git?.close()
        }
        if( !keepGit )
            FileHelper.deletePath(target.resolve('.git'))
        return head
    }

    /** Shallow-fetch the commit; servers that refuse unadvertised commits get a full clone */
    private static Git fetchCommit(String repositoryUrl, String commitId, Path target, int depth) {
        Git git = Git.init().setDirectory(target.toFile()).call()
        try {
            try {
                git.remoteAdd().setName('origin').setUri(new URIish(repositoryUrl)).call()
                git.fetch().setRemote('origin').setRefSpecs(new RefSpec(commitId)).setDepth(depth).call()
            }
            catch( GitAPIException | JGitInternalException e ) {
                log.debug "Cannot fetch commit ${commitId} of ${repositoryUrl}, falling back to a full clone: ${e.message}"
                git.close()
                FileHelper.deletePath(target)
                git = Git.cloneRepository().setURI(repositoryUrl).setDirectory(target.toFile()).call()
            }
            git.checkout().setName(commitId).call()
            return git
        }
        catch( Throwable e ) {
            git.close()
            throw e
        }
    }
}
