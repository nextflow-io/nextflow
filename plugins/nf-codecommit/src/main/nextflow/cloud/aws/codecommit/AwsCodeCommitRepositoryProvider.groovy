/*
 * Copyright 2020-2021, Seqera Labs
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
 *
 */

package nextflow.cloud.aws.codecommit

import com.amazonaws.services.codecommit.AWSCodeCommit
import com.amazonaws.services.codecommit.AWSCodeCommitClientBuilder
import com.amazonaws.services.codecommit.model.GetFileRequest
import com.amazonaws.services.codecommit.model.GetRepositoryRequest
import com.amazonaws.services.codecommit.model.RepositoryMetadata
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.scm.ProviderConfig
import nextflow.scm.RepositoryProvider
import org.eclipse.jgit.transport.CredentialsProvider
/**
 * Implements a repository provider for AWS CodeCommit
 *
 * @author W. Lee Pang <wleepang@gmail.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsCodeCommitRepositoryProvider extends RepositoryProvider {

    AwsCodeCommitRepositoryProvider(String project, ProviderConfig config) {
        assert config instanceof AwsCodeCommitProviderConfig
        this.project = project  // expect: "codecommit-<region>/<repository>"
        this.config = config
        this.region = config.region
        this.repositoryName = project.tokenize('/')[-1]
        this.client = AWSCodeCommitClientBuilder.standard().withRegion(region).build()
    }

    private String region
    private AWSCodeCommit client
    private String repositoryName


    /** {@inheritDoc} **/
    @Memoized
    @Override
    CredentialsProvider getGitCredentials() {
        return new AwsCodeCommitCredentialProvider()
    }

    private RepositoryMetadata getRepositoryMetadata() {
        final request = new GetRepositoryRequest()
                .withRepositoryName(repositoryName)

        return client
                .getRepository(request)
                .getRepositoryMetadata()
    }

    /** {@inheritDoc} **/
    // called by AssetManager
    // used to set credentials for a clone, pull, fetch, operation
    @Override
    boolean hasCredentials() {
        // set to true
        // uses AWS Credentials instead of username : password
        // see getGitCredentials()
        return true
    }

    /** {@inheritDoc} **/
    @Override
    String getName() { "CodeCommit" }

    /** {@inheritDoc} **/
    @Override
    String getEndpointUrl() {
        "https://git-codecommit.${region}.amazonaws.com/v1/repos/${repositoryName}"
    }

    /** {@inheritDoc} **/
    // not used, but the abstract method needs to be overridden
    @Override
    String getContentUrl( String path ) {
        throw new UnsupportedOperationException()
    }

    /** {@inheritDoc} **/
    // called by AssetManager
    @Override
    String getCloneUrl() { getEndpointUrl() }

    /** {@inheritDoc} **/
    // called by AssetManager
    @Override
    String getRepositoryUrl() { getEndpointUrl() }

    /** {@inheritDoc} **/
    // called by AssetManager
    // called by RepositoryProvider.readText()
    @Override
    protected byte[] readBytes( String path ) {

        final request = new GetFileRequest()
            .withRepositoryName(repositoryName)
            .withFilePath(path)
        if( revision )
            request.withCommitSpecifier(revision)

        try {
            return client
                    .getFile( request )
                    .getFileContent()?.array()
        }
        catch (Exception e) {
            log.debug "AWS CodeCommit unable to retrieve file: $path from repo: $repositoryName"
            return null
        }
    }

    /** {@inheritDoc} **/
    // called by AssetManager
    @Override
    void validateRepo() {
        try {
            getRepositoryMetadata()
        }
        catch( IOException e ) {
            throw new AbortOperationException("Cannot find ${getEndpointUrl()} -- Make sure a repository exists for it in AWS CodeCommit")
        }
    }

}
