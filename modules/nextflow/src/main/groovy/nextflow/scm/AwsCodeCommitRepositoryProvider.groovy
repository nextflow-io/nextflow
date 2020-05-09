package nextflow.scm

import nextflow.Global
import nextflow.cloud.aws.AmazonCloudDriver
import nextflow.exception.AbortOperationException

import groovy.util.logging.Slf4j
import groovy.transform.CompileStatic

import com.amazonaws.auth.AWSCredentialsProvider
import com.amazonaws.auth.profile.ProfileCredentialsProvider
import com.amazonaws.services.codecommit.AWSCodeCommitClient
import com.amazonaws.services.codecommit.AWSCodeCommitClientBuilder
import com.amazonaws.services.codecommit.model.GetRepositoryRequest
import com.amazonaws.services.codecommit.model.GetRepositoryResult
import com.amazonaws.services.codecommit.model.RepositoryMetadata
import com.amazonaws.services.codecommit.model.GetFileRequest

import org.springframework.cloud.config.server.support.AwsCodeCommitCredentialProvider

import org.eclipse.jgit.transport.CredentialsProvider

/**
 * Implements a repository provider for AWS CodeCommit
 *
 * @author W. Lee Pang <wleepang@gmail.com>
 */
@Slf4j
final class AwsCodeCommitRepositoryProvider extends RepositoryProvider {

    AwsCodeCommitRepositoryProvider(String project, ProviderConfig config=null, AWSCodeCommitClient client=null) {
        this.driver = new AmazonCloudDriver()

        this.project = project  // expect: "codecommit:[region]://<repository>"
        this.config = config ?: new ProviderConfig('codecommit')
        this.region = getRepositoryRegion()
        this.repositoryName = getRepositoryName()
        this.profile = getRepositoryProfile()
        
        this.client = client ?: getClient()

    }

    private AmazonCloudDriver driver
    private String region
    private AWSCodeCommitClient client
    private String repositoryName
    private String profile
    private RepositoryMetadata repositoryMetadata

    /** {@inheritDoc} **/
    @Override
    CredentialsProvider getGitCredentials() {
        def provider = new AwsCodeCommitCredentialProvider()
        provider.setAwsCredentialProvider( getAwsCredentialsProvider() )
        return provider
    }

    private String getRepositoryName() {
        def result = project
            .replaceAll("codecommit:.*?//", "")
            .replaceAll(".*?@", "")
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

    private String getRepositoryRegion() {
        def _region = driver.region

        // use the repository region if specified
        def result = project
            .replaceAll("codecommit:[:]*", "")
            .replaceAll("[:]*//.*", "")
        if ( result != "" ) {
            _region = result
        }

        log.debug "project region: $_region"
        return _region
    }

    private AWSCredentialsProvider getAwsCredentialsProvider() {
        if ( profile ) {
            // use locally configured profile
            return new ProfileCredentialsProvider( profile )
        }

        return driver.getCredentialsProvider0()
    }

    private AWSCodeCommitClient getClient() {

        def builder = AWSCodeCommitClientBuilder.standard()
        builder.setCredentials( getAwsCredentialsProvider() )
        if ( region ) {
            builder.setRegion(region)
        }

        return builder.build()
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
    String getContentUrl( String path ) {}

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