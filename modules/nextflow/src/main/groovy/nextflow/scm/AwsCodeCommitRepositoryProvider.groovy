package nextflow.scm

import nextflow.Global
import nextflow.cloud.aws.AmazonCloudDriver
import nextflow.exception.AbortOperationException

import groovy.util.logging.Slf4j
import groovy.transform.CompileStatic

import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.auth.credentials.StaticCredentialsProvider
import software.amazon.awssdk.auth.credentials.AwsBasicCredentials
import software.amazon.awssdk.auth.credentials.AwsSessionCredentials
import software.amazon.awssdk.services.codecommit.CodeCommitClient
import software.amazon.awssdk.services.codecommit.model.GetRepositoryRequest
import software.amazon.awssdk.services.codecommit.model.GetFileRequest
import software.amazon.awssdk.services.codecommit.model.RepositoryMetadata

/**
 * Implements a repository provider for AWS CodeCommit
 *
 * @author W. Lee Pang <wleepang@gmail.com>
 */
@Slf4j
final class AwsCodeCommitRepositoryProvider extends RepositoryProvider {

    AwsCodeCommitRepositoryProvider(String project, ProviderConfig config=null) {
        this.driver = new AmazonCloudDriver()

        this.project = project  // expect: "codecommit:[region]://<repository>"
        this.config = config ?: new ProviderConfig('codecommit')
        this.region = getRepositoryRegion()
        
        this.client = getClient()
        this.repositoryName = getRepositoryName()
        this.repositoryMetadata = getRepositoryMetadata()
    }

    private AmazonCloudDriver driver
    private Region region
    private CodeCommitClient client
    private String repositoryName
    private RepositoryMetadata repositoryMetadata

    private String getRepositoryName() {
        def name = project.replaceAll("codecommit:.*?//", "")
        log.debug "project name: $name"
        return name
    }

    private Region getRepositoryRegion() {
        def _region = driver.region

        // use the repository region if specified
        def projectRegion = project.replaceAll("codecommit:[:]*", "").replaceAll("[:]*//.*", "")
        if ( projectRegion != "" ) {
            _region = projectRegion
        }

        log.debug "project region: $_region"
        return Region.of( _region )
    }

    private CodeCommitClient getClient() {

        // aws codecommit repositories can be from different regions

        def builder = CodeCommitClient.builder()

        if ( region ) {
            builder.region(region)
        }

        if ( driver.getAccessKey() && driver.getSecretKey() ) {
            // use use aws config scope for credentials
            if ( driver.getSessionToken() ) {
                client = builder
                    .credentialsProvider(StaticCredentialsProvider.create(
                        AwsSessionCredentials.create(
                            driver.getAccessKey(),
                            driver.getSecretKey(),
                            driver.getSessionToken()
                        )
                    ))
                    .build()
            } else {
                client = builder
                    .credentialsProvider(StaticCredentialsProvider.create(
                        AwsBasicCredentials.create(
                            driver.getAccessKey(),
                            driver.getSecretKey()
                        )
                    ))
                    .build()
            }

        } else {
            // otherwise use the default credentials chain
            // this is the recommended way to provide AWS credentials
            client = builder.build()
        }

        return client
    }

    private RepositoryMetadata getRepositoryMetadata() {
        def response = client.getRepository(
            GetRepositoryRequest.builder()
                .repositoryName( repositoryName )
                .build()
        )
        return response.repositoryMetadata()
    }

    /** {@inheritDoc} **/
    // called by AssetManager
    // used to set credentials for a clone, pull, fetch, operation
    @Override
    boolean hasCredentials() {
        // set to false
        // use AWS Credentials instead of username : password
        return false
    }

    // called by AssetManager
    // RepositoryProvider setCredentials(String userName, String password)

    /** {@inheritDoc} **/
    @Override
    String getName() { "CodeCommit" }

    /** {@inheritDoc} **/
    @Override
    String getEndpointUrl() {
        //"${config.endpoint}/v1/repos/${project}"
    }

    /** {@inheritDoc} **/
    // not used, but the abstract method needs to be overriden
    @Override
    String getContentUrl( String path ) {}

    /** {@inheritDoc} **/
    // called by AssetManager
    @Override
    String getCloneUrl() { repositoryMetadata.cloneUrlHttp() }

    /** {@inheritDoc} **/
    // called by AssetManager
    @Override
    String getRepositoryUrl() { repositoryMetadata.cloneUrlHttp() }
    
    /** {@inheritDoc} **/
    // called by AssetManager
    // called by RepositoryProvider.readText()
    @Override
    protected byte[] readBytes( String path ) {
        def response = client.getFile(
            GetFileRequest.builder()
                .repositoryName( repositoryName )
                .filePath(path)
                .build()
        )

        def contents = response.fileContent().asUtf8String()
        return contents.bytes
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