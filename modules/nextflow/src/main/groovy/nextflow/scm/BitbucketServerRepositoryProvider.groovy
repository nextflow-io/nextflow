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
import nextflow.exception.AbortOperationException
/**
 * Implements a repository provider for the BitBucket service
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
final class BitbucketServerRepositoryProvider extends RepositoryProvider {
    private String repository

    BitbucketServerRepositoryProvider(String name, ProviderConfig config=null) {
        def parts = name.split('/') as List<String>
        if( parts.size()!=2 )
            throw new AbortOperationException("Not a valid project name: $name")

        def repo = parts[-1]
        def project = parts[0]
        println("Split project: ${project}, repo: ${repo}")
        this.project = project
        this.repository = repo
        this.config = config ?: new ProviderConfig('bitbucketserver')
    }

    /** {@inheritDoc} */
    @Override
    String getName() { "BitBucketServer" }

    @Override
    String getEndpointUrl() {
        println("Called BitbucketServerRepositoryProvider::getEndpointUrl: ${config.endpoint}/rest/api/1.0/projects/${project}")
        return "${config.endpoint}/rest/api/1.0/projects/${project}/repos/${repository}"
    }

    @Override
    String getContentUrl( String path ) {
        println("Called getContentUrl: ${config.endpoint}/rest/api/1.0/projects/${project}/repos/${repository}/raw/${path}")
        // "${config.endpoint}/projects/$project/src/${getMainBranch()}/$path"
        return  "${config.endpoint}/rest/api/1.0/projects/${project}/repos/${repository}/raw/${path}"    
    }

    private String getMainBranchUrl() {
        println("Called getMainBranchUrl:")
        println("${config.endpoint}/rest/api/1.0/projects/$project/repos/${repository}/branches/default")
        return  "${config.endpoint}/rest/api/1.0/projects/$project/repos/${repository}/branches/default"
    }

    String getMainBranch() {
        println("Called getMainBranch:")
        return invokeAndParseResponse(getMainBranchUrl()) ?. mainbranch ?. name
    }

    @Override
    String getCloneUrl() {
        println("called BitbucketServerRepositoryProvider::getCloneUrl")
        Map response = invokeAndParseResponse( getEndpointUrl() )
        println("BitbucketServerRepositoryProvider::getCloneUrl::response ${response}")

        println("BitbucketServerRepositoryProvider::getCloneUrl::response.scm ${response.scm}")

        if( response?.scmId != "git" ){
            throw new AbortOperationException("Bitbucket Server repository at ${getRepositoryUrl()} is not supporting Git")
        }

        def result = response?.links?.clone?.find{ it.name == "http" } as Map
        if( !result )
            throw new IllegalStateException("Missing clone URL for: $project")

        println("BitbucketRepositoryProvider::getCloneUrl::result.href ${result.href}")

        return result.href
    }

    @Override
    String getRepositoryUrl() {
        return "${config.server}/scm/${project}/${repository}"
    }

    @Override
    byte[] readBytes(String path) {
        println("called BitbucketServerRepositoryProvider::readBytes with path: ${path}")
        def url = getContentUrl(path)
        invoke(url)?.getBytes()
    }
}