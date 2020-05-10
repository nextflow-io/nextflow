/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic

/**
 * Implements a repository provider for Gitea service
 *
 * @author Akira Sekiguchi <pachiras.yokohama@gmail.com>
 */
@CompileStatic
final class GiteaRepositoryProvider extends RepositoryProvider {

    GiteaRepositoryProvider(String project, ProviderConfig config=null) {
        this.project = project
        this.config = config ?: new ProviderConfig('gitea')
    }

    /** {@inheritDoc} */
    @Override
    String getName() { "Gitea" }

    /** {@inheritDoc} */
    @Override
    String getEndpointUrl() {
        "${config.endpoint}/repos/${project}"
    }

    @Override
    @CompileDynamic
    List<BranchInfo> getBranches() {
        // https://try.gitea.io/api/swagger#/repository/repoListBranches
        final url = "${config.endpoint}/repos/${project}/branches"
        this.<BranchInfo>invokeAndResponseWithPaging(url, { Map branch -> new BranchInfo(branch.name as String, branch.commit?.id as String) })
    }

    @Override
    @CompileDynamic
    List<TagInfo> getTags() {
        // https://try.gitea.io/api/swagger#/repository/repoListTags
        final url = "${config.endpoint}/repos/${project}/tags"
        this.<TagInfo>invokeAndResponseWithPaging(url, { Map tag -> new TagInfo(tag.name as String, tag.commit?.sha as String) })
    }

    /** {@inheritDoc} */
    @Override
    String getContentUrl( String path ) {
        "${config.endpoint}/repos/$project/raw/$path"
    }

    /** {@inheritDoc} */
    @Override
    String getCloneUrl() {
        Map response = invokeAndParseResponse( getEndpointUrl() )

        def result = response.get('clone_url')
        if( !result )
            throw new IllegalStateException("Missing clone URL for: $project")

        return result
    }

    /** {@inheritDoc} */
    @Override
    String getRepositoryUrl() {
        "${config.server}/$project"
    }

    /** {@inheritDoc} */
    @Override
    byte[] readBytes(String path) {

        def url = getContentUrl(path)
        def contents = invoke(url)
        return contents?.getBytes()
    }

}
