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
import nextflow.config.Manifest
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.api.errors.RefNotFoundException
import org.eclipse.jgit.lib.Ref

/**
 * Strategy interface for different repository management approaches.
 * Implementations handle the specifics of how pipelines are stored and accessed locally.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
interface RepositoryStrategy {

    /**
     * Download or update the repository.
     * For legacy strategy, this clones the repo directly to the local path.
     * For multi-revision strategy, this creates a bare repo and shared clones.
     *
     * @param revision The revision to download (branch, tag, or commit SHA)
     * @param deep Optional depth for shallow clones
     * @return Status message describing what was done
     */
    String download(String revision, Integer deep, Manifest manifest)

    /**
     * Checkout a specific revision and returns exception if revision is not found locally
     * @param revision The revision to be checked out
     */
    void tryCheckout(String revision, Manifest manifest) throws RefNotFoundException

    /**
     * Fetch and checkout a specific revision
     * @param revision The revision to be checked out
     */
    Ref checkoutRemoteBranch(String revision, Manifest manifest)

    /**
     * @return The local path where the worktree and repository files are accessible
     */
    File getLocalPath()

    /**
     * @return The Git repository object for repository operations
     */
    Git getGit()

    /**
     * @return True if the working directory has no uncommitted changes
     */
    boolean isClean()

    /**
     * @return The default branch name from the remote repository
     */
    String getRemoteDefaultBranch()

    /**
     * @return List of all branches in the repository
     */
    List<Ref> getBranchList()

    /**
     * @return List of all tags in the repository
     */
    List<Ref> getTagList()

    /**
     * Query remote repository for refs
     *
     * @param tags If true, query for tags; otherwise query for branches
     * @return Map of ref names to Ref objects
     */
    Map<String, Ref> lsRemote(boolean tags)

    /**
     * Peel a ref (resolves annotated tags to their target)
     *
     * @param ref The ref to peel
     * @return The peeled ref
     */
    Ref peel(Ref ref)

    /**
     * @return The local git config file
     */
    File getLocalGitConfig()

    /**
     * @return The local asset repository URL
     */
    String getGitRepositoryUrl()

    /**
     * @return The current revision in the local repository
     */
    String getCurrentRevision()

    /**
     * @return The list of current downloaded commits in the local repository
     */
    List<String> listDownloadedCommits()

    /**
     * @return The current revision and Name
     */
    AssetManager.RevisionInfo getCurrentRevisionAndName()

    /**
     * Drop local copy of a repository. If revision is specified, only removes the specified revision
     * @param revision
     */
    void drop(String revision, boolean force)

    /**
     * Close any open resources (e.g., Git objects)
     */
    void close()

    /**
     * Set the localPath
     * @param file
     */
    void setLocalPath(File file)

    /**
     * Set project
     * @param projectName
     */
    void setProject(String projectName)

    /**
     * @return return the local path where project repository and revisions are download
     */
    File getProjectPath()

    /**
     * Set repository provider
     * @param provider
     */
    void setProvider(RepositoryProvider provider)
}
