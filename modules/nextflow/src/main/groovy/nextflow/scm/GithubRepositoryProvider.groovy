/*
 * Copyright 2013-2024, Seqera Labs
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
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.SysEnv

/**
 * Implements a repository provider for GitHub service
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class GithubRepositoryProvider extends RepositoryProvider {

    GithubRepositoryProvider(String project, ProviderConfig config=null) {
        this.project = project
        this.config = config ?: new ProviderConfig('github')
    }

    /** {@inheritDoc} */
    @Override
    String getName() { "GitHub" }

    /** {@inheritDoc} */
    @Override
    String getEndpointUrl() {
        return "${config.endpoint}/repos/${project}"
    }

    @Override
    boolean hasCredentials() {
        super.hasCredentials() ?: SysEnv.containsKey('GITHUB_TOKEN')
    }

    @Override
    String getUser() {
        super.getUser() ?: SysEnv.get('GITHUB_TOKEN')
    }

    @Override
    String getPassword() {
        super.getPassword() ?: (SysEnv.containsKey('GITHUB_TOKEN') ? 'x-oauth-basic' : null)
    }

    /** {@inheritDoc} */
    @Override
    String getContentUrl( String path ) {
        // see
        // https://docs.github.com/en/rest/reference/repos#get-repository-content
        //
        def result = "${config.endpoint}/repos/$project/contents/$path"
        if( revision )
            result += "?ref=$revision"
        return result
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

    @Override
    @CompileDynamic
    @Memoized
    List<BranchInfo> getBranches() {
        // https://developer.github.com/v3/repos/branches/#list-branches
        final url = "${config.endpoint}/repos/$project/branches"
        this.<BranchInfo>invokeAndResponseWithPaging(url, { Map branch -> new BranchInfo(branch.name as String, branch.commit?.sha as String) })
    }

    @Override
    @CompileDynamic
    @Memoized
    List<TagInfo> getTags() {
        // https://developer.github.com/v3/repos/#list-tags
        final url = "${config.endpoint}/repos/$project/tags"
        this.<TagInfo>invokeAndResponseWithPaging(url, { Map tag -> new TagInfo(tag.name as String, tag.commit?.sha as String)})
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
        Map response  = invokeAndParseResponse(url)
        response.get('content')?.toString()?.decodeBase64()

    }

}
