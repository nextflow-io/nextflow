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
import org.eclipse.jgit.api.Git

/**
 * Fetches the exact tree of a pinned commit from a public Git repository.
 *
 * Only public repositories are supported (no credentials). Switch to
 * AssetManager if private registries arrive.
 */
@CompileStatic
class PipelineGitFetcher {

    void fetch(String repositoryUrl, String revision, String commitId, Path target) {
        final isSha = revision ==~ /[0-9a-fA-F]{40}/
        def clone = Git.cloneRepository()
            .setURI(repositoryUrl)
            .setDirectory(target.toFile())
        if( !isSha )
            clone.setBranch(revision).setDepth(1)
        def git = clone.call()
        try {
            if( isSha )
                git.checkout().setName(commitId).call()
            final head = git.repository.resolve('HEAD')?.name()
            if( head != commitId )
                throw new AbortOperationException("Revision '$revision' of $repositoryUrl does not match the pinned commit $commitId")
        }
        finally {
            git.close()
        }
        target.resolve('.git').toFile().deleteDir()
    }
}
