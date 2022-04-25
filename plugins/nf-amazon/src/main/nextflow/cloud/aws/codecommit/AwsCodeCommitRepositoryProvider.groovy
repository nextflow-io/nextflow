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
import com.amazonaws.services.codecommit.AWSCodeCommitClient
import com.amazonaws.services.codecommit.model.GetFileRequest
import com.amazonaws.services.codecommit.model.GetRepositoryRequest
import com.amazonaws.services.codecommit.model.RepositoryMetadata
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.cloud.aws.AmazonClientFactory
import nextflow.exception.AbortOperationException
import nextflow.scm.ProviderConfig
import nextflow.scm.RepositoryProvider
import org.eclipse.jgit.transport.CredentialsProvider

import java.util.regex.Matcher

/**
 * Implements a repository provider for AWS CodeCommit
 *
 * @author W. Lee Pang <wleepang@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsCodeCommitRepositoryProvider extends RepositoryProvider {

    private AmazonClientFactory clientFactory

    AwsCodeCommitRepositoryProvider(String project, String url, ProviderConfig config=null,  AWSCodeCommitClient client=null) {

        this.project = project  // expect: "codecommit:[region]://<repository>"
        this.config = config ?: new ProviderConfig('codecommit')
        this.region = getRepositoryRegion(url)
        this.repositoryName = getRepositoryName()
        this.profile = getRepositoryProfile()
        this.clientFactory = AmazonClientFactory.instance([region:region])
        this.client = client ?: clientFactory.getCodeCommitClient()

    }

    private String region
    private AWSCodeCommit client
    private String repositoryName
    private String profile


    /** {@inheritDoc} **/
    @Override
    CredentialsProvider getGitCredentials() {
        def provider = new AwsCodeCommitCredentialProvider()
        provider.setAwsCredentialsProvider( clientFactory.getCredentialsProvider() )
        return provider
    }

    private String getRepositoryName() {
        def result = project
                .replaceAll("codecommit:.*?//", "")
                .replaceAll(".*?@", "")
                .replaceAll("repos/", "")
        log.debug "project name: $result"
        return result
    }

    private String getRepositoryProfile() {

        def result = null

        if ( project.indexOf('@') != -1 ) {
            result = project
                    .replaceAll("codecommit:.*?//", "")
                    .replaceAll("@.*", "")

            log.debug "project profile: $result"

        }

        return result
    }

    private String getRepositoryRegion(String url) {
        def region = Global.getAwsRegion()

        // use the repository region if specified
        def pattern = /codecommit(:{2}|.)(?<region>[a-z1-9-]+)(:|.)/
        def matcher = url =~ pattern
        if(matcher.find()){
           region = matcher.group("region")
        }

        log.debug "AWS CodeCommit project region: $region"
        return region
    }


    private RepositoryMetadata getRepositoryMetadata() {
        def request = new GetRepositoryRequest()
        request.setRepositoryName( repositoryName)

        def response = client.getRepository( request )
        return response.getRepositoryMetadata()
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
        "https://git-codecommit.${region}.amazonaws.com/v1/repos/${getRepositoryName()}"
    }

    /** {@inheritDoc} **/
    // not used, but the abstract method needs to be overriden
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

        def request = new GetFileRequest()
        request.setRepositoryName( repositoryName )
        request.setFilePath( path )

        def response = client.getFile( request )

        def contents = response.getFileContent().array()
        return contents
    }

    /** {@inheritDoc} **/
    // called by AssetManager
    @Override
    void validateRepo() {
        try {
            getRepositoryMetadata()
        }
        catch( IOException e ) {
            throw new AbortOperationException("Cannot find `$project` -- Make sure a repository exists for it in AWS CodeCommit")
        }
    }

}
