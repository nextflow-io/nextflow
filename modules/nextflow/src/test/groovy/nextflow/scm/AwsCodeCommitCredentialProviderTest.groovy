package nextflow.scm

import spock.lang.Requires
import spock.lang.Specification

import nextflow.Global

import org.eclipse.jgit.transport.URIish

/**
 *
 * @author W. Lee Pang <wleepang@gmail.com>
 */
class AwsCodeCommitCredentialProviderTest extends Specification {
    def provider
    def uri
    def awsSecretKey

    def "should return sig-v4 password" () {
        
        when:
        provider = new AwsCodeCommitCredentialProvider()
        uri = new URIish("https://git-codecommit.us-west-2.amazonaws.com/v1/repos/project")
        awsSecretKey = "abc123xyz"
        
        then:
        provider.calculateCodeCommitPassword(uri, awsSecretKey, new Date(0)) == "19700101T000000Z105088f0165351288fd9bad861105f60ce80fc5da9d824b68e0601155157c657"

    }
}