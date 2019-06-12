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