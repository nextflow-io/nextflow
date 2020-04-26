package nextflow.scm

import groovy.util.logging.Slf4j
import groovy.transform.CompileStatic

import software.amazon.awssdk.services.codecommit.CodeCommitClient

/**
 * Implements a repository provider for AWS CodeCommit
 *
 * @author W. Lee Pang <wleepang@gmail.com>
 */
@Slf4j
@CompileStatic
final class AwsCodeCommitRepositoryProvider extends RepositoryProvider {

    AwsCodeCommitRepositoryProvider(String project, ProviderConfig config=null) {
        this.project = project
        this.config = config ?: new ProviderConfig('codecommit')
    }

    /** {@inheritDoc} **/
    @Override
    boolean hasCredentials() {
        // @TODO
        // use default credentials chain first
        // last resort, use keys in scm config file
    }

    /** {@inheritDoc} **/
    @Override
    String getName() { "AWS CodeCommit" }

    /** {@inheritDoc} **/
    @Override
    String getEndpointUrl() {
        "${config.endpoint}/v1/repos/${project}"
    }

    /** not implemented
     * String getContentUrl() {}
     */

    /** {@inheritDoc} **/
    @Override
    String getCloneUrl() {
        // @TODO
    }

    /** {@inheritDoc} **/
    @Override
    String getRepositoryUrl() {}
    
    /** {@inheritDoc} **/
    @Override
    protected String invoke( String api ) {}

    /** {@inheritDoc} **/
    @Override
    protected byte[] readBytes( String path ) {
        // called by readText()
    }

    // /** {@inheritDoc} **/
    // @Override
    // String readText( String path ) {}

    // /** {@inheritDoc} **/
    // @Override
    // void validateFor() {}

    // /** {@inheritDoc} **/
    // @Override
    // void validateRepo() {}

}