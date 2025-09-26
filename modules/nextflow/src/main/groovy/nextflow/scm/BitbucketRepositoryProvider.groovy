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

import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import org.eclipse.jgit.transport.CredentialsProvider
import org.eclipse.jgit.transport.UsernamePasswordCredentialsProvider

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

    @Override
    protected String[] getAuth() {
        if (!hasCredentials()) {
            return null
        }

        String secret = getToken() ?: getPassword()
        String authString = "${getUser()}:${secret}".bytes.encodeBase64().toString()

        return new String[] { "Authorization", "Basic " + authString }
    }

    @Override
    boolean hasCredentials() {
        return getUser() && (getPassword() || getToken())
    }

    /** {@inheritDoc} */
    @Override
    CredentialsProvider getGitCredentials() {
        return new UsernamePasswordCredentialsProvider(getUser(), getToken() ?: getPassword())
    }

    /** {@inheritDoc} */
    @Override
    String getName() { "BitBucket" }

    @Override
    String getEndpointUrl() {
        "${config.endpoint}/2.0/repositories/${project}"
    }

    @Override
    String getContentUrl( String path ) {
        final ref = revision ? getRefForRevision(revision) : getMainBranch()
        return "${config.endpoint}/2.0/repositories/$project/src/$ref/$path"
    }

    private String getMainBranchUrl() {
        "${config.endpoint}/2.0/repositories/$project"
    }

    String getMainBranch() {
        invokeAndParseResponse(getMainBranchUrl()) ?. mainbranch ?. name
    }

    private String getRefForRevision(String revision){
        // when the revision does not contain a slash just use as it is
        if( !revision.contains('/') ) {
            return revision
        }
        // when the revision contains a slash, it's needed to find the corresponding branch or tag name
        // see
        // https://github.com/nextflow-io/nextflow/issues/4599
        // https://jira.atlassian.com/browse/BCLOUD-20223
        try {
            return getRefForRevision0(revision, 'branches')
        }
        catch (Exception e1) {
            // if an error is reported it may be a tag name instead of a branch
            try {
                return getRefForRevision0(revision, 'tags')
            }
            // still failing, raise the previous error
            catch (Exception e2) {
                throw e1
            }
        }
    }

    private String getRefForRevision0(String revision, String type){
        final resp = invokeAndParseResponse("${config.endpoint}/2.0/repositories/$project/refs/$type/$revision")
        return resp?.target?.hash
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
        final url = "$config.endpoint/2.0/repositories/$project/refs/tags"
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
        final url = "$config.endpoint/2.0/repositories/$project/refs/branches"
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

        /**
         * The clone URL when an API token is used is of the form: https://{bitbucket_username}:{api_token}@bitbucket.org/{workspace}/{repository}.git
         * @see <a href="https://support.atlassian.com/bitbucket-cloud/docs/using-api-tokens/">Using API tokens</a>
         */
        String cloneUrl = result.href
        return getToken() ? cloneUrl.replace('@', ":${getToken()}@") : cloneUrl
    }

    @Override
    String getRepositoryUrl() {
        return "${config.server}/$project"
    }

    @Override
    byte[] readBytes(String path) {
        final url = getContentUrl(path)
        return invokeBytes(url)
    }

    /** {@inheritDoc} */
    @Override
    List<RepositoryEntry> listDirectory(String path, int depth) {
        throw new UnsupportedOperationException("Directory listing not yet implemented for BitBucket")
    }
}
