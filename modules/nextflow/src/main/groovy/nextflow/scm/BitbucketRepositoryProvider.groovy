/*
 * Copyright 2020, Seqera Labs
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

import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
/**
 * Implements a repository provider for the BitBucket service
 *
 * Nextflow uses Basic authentication for API access. It requires
 * the use of App password in place of the user password.
 *
 * See more:
 *   https://support.atlassian.com/bitbucket-cloud/docs/app-passwords/
 *   https://bitbucket.org/account/settings/app-passwords/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
final class BitbucketRepositoryProvider extends RepositoryProvider {

    BitbucketRepositoryProvider(String project, ProviderConfig config=null) {
        this.project = project
        this.config = config ?: new ProviderConfig('bitbucket')
    }

    /** {@inheritDoc} */
    @Override
    String getName() { "BitBucket" }

    @Override
    String getEndpointUrl() {
        "${config.endpoint}/api/2.0/repositories/${project}"
    }

    @Override
    String getContentUrl( String path ) {
        "${config.endpoint}/api/2.0/repositories/$project/src/${getMainBranch()}/$path"
    }

    private String getMainBranchUrl() {
        "${config.endpoint}/api/2.0/repositories/$project"
    }

    String getMainBranch() {
        invokeAndParseResponse(getMainBranchUrl()) ?. mainbranch ?. name
    }

    @Memoized
    protected <T> List<T> invokeAndResponseWithPaging(String url, Closure<T> parse) {
        final result = new ArrayList<T>(50)
        do {
            final resp = invokeAndParseResponse(url)
            final tags = resp.values as List<Map>
            for( Map entry : tags ) {
                parse.call(entry)
            }
            // next iteration
            // https://developer.atlassian.com/bitbucket/api/2/reference/meta/pagination
            url = resp.next
        }
        while( url != null )
        return result
    }

    /**
     * Fetch the repository tags
     * https://developer.atlassian.com/bitbucket/api/2/reference/resource/repositories/%7Bworkspace%7D/%7Brepo_slug%7D/refs/tags
     *
     * @return A list of {@link TagInfo}
     */
    @Override
    List<TagInfo> getTags() {
        final result = new ArrayList<TagInfo>()
        final url = "$config.endpoint/api/2.0/repositories/$project/refs/tags"
        final mapper = { Map entry -> result.add( new TagInfo(entry.name, entry.target?.hash) ) }
        invokeAndResponseWithPaging(url, mapper)
        return result
    }

    /**
     * https://developer.atlassian.com/bitbucket/api/2/reference/resource/repositories/%7Bworkspace%7D/%7Brepo_slug%7D/refs/branches
     *
     * @return A list of {@link BranchInfo}
     */
    @Override
    List<BranchInfo> getBranches() {
        final result = new ArrayList<BranchInfo>()
        final url = "$config.endpoint/api/2.0/repositories/$project/refs/branches"
        final mapper = { Map entry -> result.add( new BranchInfo(entry.name, entry.target?.hash) ) }
        invokeAndResponseWithPaging(url, mapper)
        return result
    }

    @Override
    String getCloneUrl() {
        Map response = invokeAndParseResponse( getEndpointUrl() )

        if( response?.scm != "git" ){
            throw new AbortOperationException("Bitbucket repository at ${getRepositoryUrl()} is not supporting Git")
        }

        def result = response?.links?.clone?.find{ it.name == "https" } as Map
        if( !result )
            throw new IllegalStateException("Missing clone URL for: $project")

        return result.href
    }

    @Override
    String getRepositoryUrl() {
        return "${config.server}/$project"
    }

    @Override
    byte[] readBytes(String path) {

        def url = getContentUrl(path)
        invoke(url)?.getBytes()
    }
}
