package nextflow.scm

import spock.lang.Requires
import spock.lang.Specification

import nextflow.Global

import com.amazonaws.services.codecommit.AWSCodeCommitClient
import com.amazonaws.auth.profile.ProfileCredentialsProvider

/**
 *
 * @author W. Lee Pang <wleepang@gmail.com>
 */
class AwsCodeCommitRepositoryProviderTest extends Specification {
    def provider
    
    ProviderConfig config = new ProviderConfig('codecommit')
    
    def setup() {
        GroovySpy(Global, global: true)
        _ * Global.getAwsCredentials(*_) >> ['abc', 'xyz']
        _ * Global.getAwsRegion(*_) >> 'us-west-2'
    }

    def "should return project regionalized endpoint url" () {

        when:
        provider = new AwsCodeCommitRepositoryProvider('codecommit::us-east-1://project', config)

        then:
        provider.getEndpointUrl() == "https://git-codecommit.us-east-1.amazonaws.com/v1/repos/project"
    
    }
    
    def "should return default regionalized endpoint url" () {
        when:
        provider = new AwsCodeCommitRepositoryProvider('codecommit://project', config)

        then:
        provider.getEndpointUrl() == "https://git-codecommit.us-west-2.amazonaws.com/v1/repos/project"
    }

    def "should ignore profile in project name" () {
        when:
        provider = new AwsCodeCommitRepositoryProvider('codecommit://profile@project', config)

        then:
        provider.getEndpointUrl() == "https://git-codecommit.us-west-2.amazonaws.com/v1/repos/project"
    }

    def "should use profile credential provider if profile specified" () {
        when:
        GroovySpy(ProfileCredentialsProvider, global: true)
        provider = new AwsCodeCommitRepositoryProvider('codecommit://profile@project', config)

        then:
        (1.._) * new ProfileCredentialsProvider("profile")
    }

}
