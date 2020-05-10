/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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


import groovy.util.logging.Slf4j
/**
 * Implements a repository provider for GitHub service
 *
 * See https://gitlab.com/
 * https://docs.gitlab.com/ee/api/repositories.html
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class GitlabRepositoryProvider extends RepositoryProvider {

    GitlabRepositoryProvider(String project, ProviderConfig config=null) {
        this.project = project
        this.config = config ?: new ProviderConfig('gitlab')
    }

    final String getProjectName() {
        URLEncoder.encode(project,'utf-8')
    }

    @Override
    protected void auth( URLConnection connection ) {
        if( config.token ) {
            // set the token in the request header
            connection.setRequestProperty("PRIVATE-TOKEN", config.token)
        }
    }

    @Override
    String getName() { "GitLab" }

    @Override
    String getEndpointUrl() {
        return "${config.endpoint}/api/v4/projects/${getProjectName()}"
    }

    String getDefaultBranch() {
        def result = invokeAndParseResponse(getEndpointUrl()) ?. default_branch
        if( !result ) {
            log.debug "Unable to fetch repo default branch. Using `master` branch -- See https://gitlab.com/gitlab-com/support-forum/issues/1655#note_26132691"
            return 'master'
        }
        return result
    }

    @Override
    List<BranchInfo> getBranches() {
        // https://docs.gitlab.com/ee/api/branches.html
        final url = "${config.endpoint}/api/v4/projects/${getProjectName()}/repository/branches"
        this.<BranchInfo>invokeAndResponseWithPaging(url, { Map branch -> new BranchInfo(branch.name as String, branch.commit?.id as String) })
    }

    @Override
    List<TagInfo> getTags() {
        // https://docs.gitlab.com/ee/api/tags.html
        final url = "${config.endpoint}/api/v4/projects/${getProjectName()}/repository/tags"
        this.<TagInfo>invokeAndResponseWithPaging(url, { Map tag -> new TagInfo(tag.name as String, tag.commit?.id as String) })
    }

    /** {@inheritDoc} */
    @Override
    String getContentUrl( String path ) {
        "${config.endpoint}/api/v4/projects/${getProjectName()}/repository/files/${path}?ref=${getDefaultBranch()}"
    }

    /** {@inheritDoc} */
    @Override
    String getCloneUrl() {
        Map response = invokeAndParseResponse( getEndpointUrl() )

        def result = response.get('http_url_to_repo')
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
        Map response  = invokeAndParseResponse(url)
        response.get('content')?.toString()?.decodeBase64()

    }
}
