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
import nextflow.Const
import nextflow.config.Manifest
import org.eclipse.jgit.lib.Constants
import org.eclipse.jgit.lib.ObjectId
import org.eclipse.jgit.lib.Ref
import org.eclipse.jgit.lib.Repository

/**
 * Abstract base class providing common repository operations shared by all strategies.
 * Implements the Template Method pattern where concrete strategies override specific
 * behaviors while sharing common helper methods.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
abstract class AbstractRepositoryStrategy implements RepositoryStrategy {

    protected static final String REMOTE_REFS_ROOT = "refs/remotes/origin/"
    protected static final String REMOTE_DEFAULT_HEAD = REMOTE_REFS_ROOT + "HEAD"

    /**
     * Context providing access to shared resources and configuration
     */
    protected final RepositoryContext context

    AbstractRepositoryStrategy(RepositoryContext context) {
        this.context = context
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
        }
        else {
            return new AssetManager.RevisionInfo(head.objectId.name(), null, AssetManager.RevisionInfo.Type.COMMIT)
        }
    }

    /**
     * Find the HEAD reference
     */
    protected Ref findHeadRef() {
        getGit()?.getRepository()?.findRef(Constants.HEAD)
    }

    /**
     * Try to resolve an object id to a tag name
     */
    protected String resolveTagNameByObjectId(ObjectId objectId) {
        Collection<Ref> tags = getGit()?.getRepository()?.getRefDatabase()?.getRefsByPrefix(Constants.R_TAGS)
        return tags?.find { it.objectId == objectId }?.name
    }


    /**
     * Get the git repository URL
     */
    protected String getGitRepositoryUrl() {
        return context.provider.getCloneUrl()
    }

    /**
     * Get the default branch name
     */
    protected String getDefaultBranch(Manifest manifest) {
        return manifest?.getDefaultBranch()
            ?: getRemoteDefaultBranch()
            ?: Const.DEFAULT_BRANCH
    }


}
