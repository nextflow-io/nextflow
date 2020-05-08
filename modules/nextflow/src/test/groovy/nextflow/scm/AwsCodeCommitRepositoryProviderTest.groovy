package nextflow.scm

import spock.lang.Requires
import spock.lang.Specification

import nextflow.cloud.aws.AmazonCloudDriver

import software.amazon.awssdk.services.codecommit.CodeCommitClient

/**
 *
 * @author W. Lee Pang <wleepang@gmail.com>
 */
class AwsCodeCommitRepositoryProviderTest extends Specification {
    ProviderConfig config = new ProviderConfig('codecommit')
    def driver = Mock(AmazonCloudDriver)
    def client = Mock(CodeCommitClient)

    def setup() {
        driver.region >> "us-west-2"
    }

    def "should return project regionalized endpoint url" () {

        given:
        def provider = new AwsCodeCommitRepositoryProvider('codecommit::us-east-1://project', config, client)

        when:
        def endpoint = provider.getEndpointUrl()

        then:
        endpoint == "https://git-codecommit.us-east-1.amazonaws.com/v1/repos/project"
    
    }
    
    def "should return default regionalized endpoint url" () {
        given:
        def provider = new AwsCodeCommitRepositoryProvider('codecommit://project', config, client)

        when:
        def endpoint = provider.getEndpointUrl()

        then:
        endpoint == "https://git-codecommit.us-west-2.amazonaws.com/v1/repos/project"
    }

}
