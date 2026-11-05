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

package nextflow.pipeline

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.api.errors.GitAPIException
import org.eclipse.jgit.api.errors.JGitInternalException
import org.eclipse.jgit.transport.RefSpec

/**
 * Fetches the exact tree of a pinned commit from a public Git repository.
 *
 * Only public repositories are supported (no credentials). Switch to
 * AssetManager if private registries arrive.
 */
@CompileStatic
class PipelineGitFetcher {

    void fetch(String repositoryUrl, String revision, String commitId, Path target) {
        // a release may have no git revision: the pinned commit is then the only reference
        final byCommit = !revision || revision.equalsIgnoreCase(commitId)
        Git git = null
        try {
            git = byCommit
                ? fetchCommit(repositoryUrl, commitId, target)
                : Git.cloneRepository().setURI(repositoryUrl).setDirectory(target.toFile()).setBranch(revision).setDepth(1).call()
            final head = git.repository.resolve('HEAD')?.name()
            if( !commitId.equalsIgnoreCase(head) )
                throw new AbortOperationException("Revision '$revision' of $repositoryUrl does not match the pinned commit $commitId")
        }
        catch( GitAPIException | JGitInternalException e ) {
            throw new AbortOperationException("Cannot fetch revision '${revision ?: commitId}' of $repositoryUrl: ${e.message}", e)
        }
        finally {
            git?.close()
        }
        FileHelper.deletePath(target.resolve('.git'))
    }

    /** Shallow-fetch the commit; servers that refuse unadvertised commits get a full clone */
    private static Git fetchCommit(String repositoryUrl, String commitId, Path target) {
        Git git = Git.init().setDirectory(target.toFile()).call()
        try {
            git.fetch().setRemote(repositoryUrl).setRefSpecs(new RefSpec(commitId)).setDepth(1).call()
        }
        catch( GitAPIException | JGitInternalException e ) {
            git.close()
            FileHelper.deletePath(target)
            git = Git.cloneRepository().setURI(repositoryUrl).setDirectory(target.toFile()).call()
        }
        git.checkout().setName(commitId).call()
        return git
    }
}
