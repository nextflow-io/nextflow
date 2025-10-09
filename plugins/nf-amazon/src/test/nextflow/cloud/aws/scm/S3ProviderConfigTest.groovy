/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.cloud.aws.scm

import nextflow.Global
import nextflow.Session
import software.amazon.awssdk.auth.credentials.AwsBasicCredentials
import software.amazon.awssdk.auth.credentials.StaticCredentialsProvider
import software.amazon.awssdk.regions.Region
import spock.lang.Specification

/**
 * Tests for S3ProviderConfig
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class S3ProviderConfigTest extends Specification {

    def cleanup() {
        Global.session = null
    }

    def 'should create S3 provider config with name only'() {
        when:
        def config = new S3ProviderConfig('my-bucket')

        then:
        config.name == 'my-bucket'
        config.platform == 's3'
        config.server == 's3://my-bucket'
        config.region == Region.US_EAST_1
        config.awsCredentialsProvider != null
    }

    def 'should create S3 provider config with values map'() {
        given:
        def values = [
            platform: 's3',
            server: 's3://test-bucket',
            region: 'us-west-2'
        ]

        when:
        def config = new S3ProviderConfig('test-bucket', values)

        then:
        config.name == 'test-bucket'
        config.platform == 's3'
        config.server == 's3://test-bucket'
        config.region == Region.US_WEST_2
    }

    def 'should set region from config map'() {
        given:
        def values = [
            platform: 's3',
            region: 'eu-west-1'
        ]

        when:
        def config = new S3ProviderConfig('my-bucket', values)

        then:
        config.region == Region.EU_WEST_1
    }

    def 'should set AWS credentials from config map'() {
        given:
        def values = [
            platform: 's3',
            accessKey: 'AKIAIOSFODNN7EXAMPLE',
            secretKey: 'wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY'
        ]

        when:
        def config = new S3ProviderConfig('my-bucket', values)

        then:
        config.awsCredentialsProvider instanceof StaticCredentialsProvider
        def credentials = config.awsCredentialsProvider.resolveCredentials()
        credentials.accessKeyId() == 'AKIAIOSFODNN7EXAMPLE'
        credentials.secretAccessKey() == 'wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY'
    }

    def 'should set defaults from Global session AWS config'() {
        given:
        def sessionConfig = [
            aws: [
                region: 'ap-southeast-1',
                accessKey: 'ASIA123456789EXAMPLE',
                secretKey: 'testSecretKey123'
            ]
        ]
        Global.session = Mock(Session) {
            getConfig() >> sessionConfig
        }

        when:
        def config = new S3ProviderConfig('my-bucket')

        then:
        config.region == Region.AP_SOUTHEAST_1
        config.awsCredentialsProvider instanceof StaticCredentialsProvider
    }

    def 'should override session config with provider values'() {
        given:
        def sessionConfig = [
            aws: [
                region: 'us-east-1',
                accessKey: 'ASIA111111111EXAMPLE',
                secretKey: 'sessionSecret'
            ]
        ]
        Global.session = Mock(Session) {
            getConfig() >> sessionConfig
        }
        def values = [
            platform: 's3',
            region: 'eu-central-1',
            accessKey: 'AKIAIOSFODNN7EXAMPLE',
            secretKey: 'providerSecret'
        ]

        when:
        def config = new S3ProviderConfig('my-bucket', values)

        then:
        config.region == Region.EU_CENTRAL_1
        def credentials = config.awsCredentialsProvider.resolveCredentials()
        credentials.accessKeyId() == 'AKIAIOSFODNN7EXAMPLE'
        credentials.secretAccessKey() == 'providerSecret'
    }

    def 'should resolve project name correctly'() {
        given:
        def config = new S3ProviderConfig('my-bucket')

        when:
        def result = config.resolveProjectName('path/to/project')

        then:
        result == 'my-bucket/path/to/project'
    }

    def 'should handle project name with nested paths'() {
        given:
        def config = new S3ProviderConfig('test-bucket')

        when:
        def result = config.resolveProjectName('org/team/project')

        then:
        result == 'test-bucket/org/team/project'
    }

    def 'should use default credentials provider when no credentials specified'() {
        when:
        def config = new S3ProviderConfig('my-bucket')

        then:
        config.awsCredentialsProvider != null
        // DefaultCredentialsProvider is used by default
    }

    def 'should handle different AWS regions'() {
        expect:
        new S3ProviderConfig('bucket', [platform: 's3', region: REGION]).region == EXPECTED

        where:
        REGION              | EXPECTED
        'us-east-1'         | Region.US_EAST_1
        'us-west-2'         | Region.US_WEST_2
        'eu-west-1'         | Region.EU_WEST_1
        'ap-northeast-1'    | Region.AP_NORTHEAST_1
        'sa-east-1'         | Region.SA_EAST_1
    }

    def 'should handle accessKey without secretKey'() {
        given:
        def values = [
            platform: 's3',
            accessKey: 'AKIAIOSFODNN7EXAMPLE'
        ]

        when:
        def config = new S3ProviderConfig('my-bucket', values)

        then:
        // Should not set static credentials provider if secretKey is missing
        config.awsCredentialsProvider != null
    }

    def 'should handle secretKey without accessKey'() {
        given:
        def values = [
            platform: 's3',
            secretKey: 'wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY'
        ]

        when:
        def config = new S3ProviderConfig('my-bucket', values)

        then:
        // Should not set static credentials provider if accessKey is missing
        config.awsCredentialsProvider != null
    }

    def 'should handle empty config map'() {
        when:
        def config = new S3ProviderConfig('my-bucket')

        then:
        config.name == 'my-bucket'
        config.platform == 's3'
        config.server == 's3://my-bucket'
        config.region == Region.US_EAST_1
    }

    def 'should handle null Global session'() {
        given:
        Global.session = null

        when:
        def config = new S3ProviderConfig('my-bucket')

        then:
        noExceptionThrown()
        config.region == Region.US_EAST_1
    }
}
