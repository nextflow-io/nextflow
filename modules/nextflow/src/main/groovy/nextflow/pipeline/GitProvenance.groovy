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
import nextflow.scm.GitUrl
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.lib.Constants
import org.eclipse.jgit.lib.Repository
import org.eclipse.jgit.revwalk.RevWalk

/**
 * Git provenance of a pipeline checkout: where it comes from and whether it is pristine
 */
@CompileStatic
class GitProvenance {

    String repositoryUrl
    String revision
    String commitId
    boolean clean

    /** @return the provenance of the checkout holding {@code dir}, or null if it is not a Git work tree */
    static GitProvenance of(Path dir) {
        final repo = PipelineBundle.openRepository(dir)
        if( repo == null )
            return null
        try {
            final head = repo.resolve(Constants.HEAD)
            if( head == null )
                return null
            final result = new GitProvenance()
            result.commitId = head.name()
            result.revision = tagAt(repo, head.name()) ?: repo.branch
            result.repositoryUrl = normalize(repo.config.getString('remote', 'origin', 'url'))
            result.clean = hasNoChanges(repo, dir)
            return result
        }
        finally {
            repo.close()
        }
    }

    private static String tagAt(Repository repo, String commitId) {
        final names = new ArrayList<String>()
        try( final walk = new RevWalk(repo) ) {
            for( final ref : repo.refDatabase.getRefsByPrefix(Constants.R_TAGS) ) {
                final id = walk.peel(walk.parseAny(ref.objectId))
                if( id.name() == commitId )
                    names.add(ref.name.substring(Constants.R_TAGS.length()))
            }
        }
        // ponytail: plain string order, not version order (v1.10 < v1.9)
        return names ? Collections.max(names) : null
    }

    /** Only changes under {@code dir} count, so unrelated edits in the repo do not block publishing */
    private static boolean hasNoChanges(Repository repo, Path dir) {
        final prefix = repo.workTree.toPath().toRealPath().relativize(dir.toRealPath()).toString()
        try( final git = new Git(repo) ) {
            final cmd = git.status()
            if( prefix )
                cmd.addPath(prefix.replace(File.separator, '/'))
            return cmd.call().isClean()
        }
    }

    /** Public https URL of the origin, or null when it cannot be expressed as one (e.g. a local path) */
    private static String normalize(String url) {
        if( !url )
            return null
        try {
            final u = new GitUrl(url)
            final path = u.path?.replaceAll('/+$', '')
            if( !(u.protocol in ['ssh', 'git', 'http', 'https']) || !u.domain || !path )
                return null
            final host = u.domain.replaceFirst(/:\d+$/, '')
            final scheme = u.protocol == 'http' ? 'http' : 'https'
            return "${scheme}://${host}/${path}".toString()
        }
        catch( IllegalArgumentException ignored ) {
            return null
        }
    }
}
