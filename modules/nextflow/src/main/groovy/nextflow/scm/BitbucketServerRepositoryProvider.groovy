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


import groovy.json.JsonSlurper
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
/**
 * Implements a repository provider for the private hosted BitBucket Server service
 *
 * See for details
 *  https://confluence.atlassian.com/bitbucketserver
 *  
 * @author Piotr Faba <piotr.faba@ardigen.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
final class BitbucketServerRepositoryProvider extends RepositoryProvider {
    private String repository
    private String givenName

    BitbucketServerRepositoryProvider(String name, ProviderConfig config=null) {
        this.givenName = name

        final parts = name.tokenize('/')
        if( parts.size() == 2 ) {
            this.project = parts[0]
            this.repository = parts[1]
        }
        else if( parts.size()==3 && parts[0]=='scm' ) {
            this.project = parts[1]
            this.repository = parts[2]
        }
        else if( parts.size() == 4 && parts[0]=='projects' && parts[2]=='repos') {
            this.project = parts[1]
            this.repository = parts[3]
        }
        else
            throw new AbortOperationException("Not a valid project name: $name")

        this.config = config ?: new ProviderConfig('bitbucketserver')
    }

    /** {@inheritDoc} */
    @Override
    String getName() { "BitBucketServer" }

    @Override
    String getEndpointUrl() {
        return "${config.endpoint}/rest/api/1.0/projects/${project}/repos/${repository}"
    }

    @Override
    String getContentUrl( String path ) {
        // see
        // https://docs.atlassian.com/bitbucket-server/rest/7.10.0/bitbucket-rest.html#idp358
        //
        def result = "${config.endpoint}/rest/api/1.0/projects/${project}/repos/${repository}/raw/${path}"
        if( revision )
            result += "?at=$revision"
        return result
    }

    private String getMainBranchUrl() {
        return  "${config.endpoint}/rest/api/1.0/projects/${project}/repos/${repository}/branches/default"
    }

    String getMainBranch() {
        return invokeAndParseResponse(getMainBranchUrl()) ?. mainbranch ?. name
    }

    @Override
    String getCloneUrl() {
        Map response = invokeAndParseResponse( getEndpointUrl() )

        if( response?.scmId != "git" ){
            throw new AbortOperationException("Bitbucket Server repository at ${getRepositoryUrl()} is not supporting Git")
        }

        def result = response?.links?.clone?.find{ it.name == "http" } as Map
        if( !result )
            throw new IllegalStateException("Missing clone URL for: ${project}")

        return result.href
    }

    @Override
    String getRepositoryUrl() {
        return "${config.server}/${givenName.stripStart('/')}"
    }

    @Override
    byte[] readBytes(String path) {
        def url = getContentUrl(path)
        invoke(url)?.getBytes()
    }

    @Override
    List<TagInfo> getTags() {
        final result = new ArrayList<TagInfo>()
        final url = "${getEndpointUrl()}/tags"
        final mapper = { Map entry -> result.add( new TagInfo(entry.displayId, entry.latestCommit) ) }
        invokeAndResponseWithPaging(url, mapper)
        return result
    }

    @Override
    List<BranchInfo> getBranches() {
        final result = new ArrayList<BranchInfo>()
        final url = "${getEndpointUrl()}/branches"
        final mapper = { Map entry -> result.add( new BranchInfo(entry.displayId, entry.latestCommit) ) }
        invokeAndResponseWithPaging(url, mapper)

        return result
    }

    @Memoized
    protected <T> List<T> invokeAndResponseWithPaging(String request, Closure<T> parse) {
        int start = 0
        final limit = 100
        final result = new ArrayList()
        while( true ) {
            final url = request + "?start=${start}&limit=${limit}"
            final response = invoke(url)
            final resp = (Map) new JsonSlurper().parseText(response)
            final list = resp.values as List
            if( !list )
                break

            for( def item : list ) {
                result.add( parse(item) )
            }

            final next = resp.nextPageStart as Integer
            if( !next || next <= start )
                break
            else
                start = next
        }
        return result
    }

}
