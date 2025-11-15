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

import nextflow.cloud.aws.scm.jgit.S3GitCredentialsProvider
import nextflow.scm.ProviderConfig
import software.amazon.awssdk.auth.credentials.DefaultCredentialsProvider
import software.amazon.awssdk.auth.credentials.StaticCredentialsProvider
import software.amazon.awssdk.regions.Region
import spock.lang.Specification

/**
 * Tests for S3RepositoryProvider
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class S3RepositoryProviderTest extends Specification {

    def 'should create S3 repository provider'() {
        given:
        def config = new S3ProviderConfig('test-bucket')

        when:
        def provider = new S3RepositoryProvider('test-bucket/project', config)

        then:
        provider.project == 'test-bucket/project'
        provider.config == config
    }

    def 'should assert config is S3ProviderConfig'() {
        given:
        def config = Mock(ProviderConfig)

        when:
        new S3RepositoryProvider('test-project', config)

        then:
        thrown(AssertionError)
    }

    def 'should return correct name'() {
        given:
        def config = new S3ProviderConfig('test-bucket')
        def provider = new S3RepositoryProvider('test-bucket/project', config)

        when:
        def name = provider.getName()

        then:
        name == 'test-bucket/project'
    }

    def 'should return correct endpoint URL'() {
        given:
        def config = new S3ProviderConfig('test-bucket')
        def provider = new S3RepositoryProvider('test-bucket/project', config)

        when:
        def url = provider.getEndpointUrl()

        then:
        url == 's3://test-bucket/project'
    }

    def 'should return correct clone URL'() {
        given:
        def config = new S3ProviderConfig('test-bucket')
        def provider = new S3RepositoryProvider('test-bucket/project', config)

        when:
        def url = provider.getCloneUrl()

        then:
        url == 's3://test-bucket/project'
    }

    def 'should return correct repository URL'() {
        given:
        def config = new S3ProviderConfig('test-bucket')
        def provider = new S3RepositoryProvider('test-bucket/project', config)

        when:
        def url = provider.getRepositoryUrl()

        then:
        url == 's3://test-bucket/project'
    }

    def 'should indicate has credentials'() {
        given:
        def config = new S3ProviderConfig('test-bucket')
        def provider = new S3RepositoryProvider('test-bucket/project', config)

        when:
        def hasCredentials = provider.hasCredentials()

        then:
        hasCredentials
    }

    def 'should throw UnsupportedOperationException for getContentUrl'() {
        given:
        def config = new S3ProviderConfig('test-bucket')
        def provider = new S3RepositoryProvider('test-bucket/project', config)

        when:
        provider.getContentUrl('path/to/file')

        then:
        thrown(UnsupportedOperationException)
    }

    def 'should get Git credentials with region'() {
        given:
        def config = new S3ProviderConfig('test-bucket', [platform: 's3', region: 'eu-west-1'])
        def provider = new S3RepositoryProvider('test-bucket/project', config)

        when:
        def credentials = provider.getGitCredentials()
        then:
        credentials instanceof S3GitCredentialsProvider
        def awsCredentials = credentials as S3GitCredentialsProvider
        awsCredentials.region == Region.EU_WEST_1
        awsCredentials.awsCredentialsProvider instanceof DefaultCredentialsProvider
    }

    def 'should get Git credentials with AWS credentials provider'() {
        given:
        def config = new S3ProviderConfig('test-bucket', [
            platform: 's3',
            region: 'us-west-2',
            accessKey: 'AKIAIOSFODNN7EXAMPLE',
            secretKey: 'wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY'
        ])
        def provider = new S3RepositoryProvider('test-bucket/project', config)

        when:
        def credentials = provider.getGitCredentials()

        then:
        credentials instanceof S3GitCredentialsProvider
        def awsCredentials = credentials as S3GitCredentialsProvider
        awsCredentials.region == Region.US_WEST_2
        awsCredentials.awsCredentialsProvider instanceof StaticCredentialsProvider
    }

    def 'should memoize Git credentials'() {
        given:
        def config = new S3ProviderConfig('test-bucket')
        def provider = new S3RepositoryProvider('test-bucket/project', config)

        when:
        def credentials1 = provider.getGitCredentials()
        def credentials2 = provider.getGitCredentials()

        then:
        credentials1.is(credentials2) // Should return same instance due to @Memoized
    }

    def 'should validate repo without throwing exception'() {
        given:
        def config = new S3ProviderConfig('test-bucket')
        def provider = new S3RepositoryProvider('test-bucket/project', config)

        when:
        provider.validateRepo()

        then:
        noExceptionThrown()
    }

    def 'should handle different project names'() {
        given:
        def config = new S3ProviderConfig('my-bucket')

        expect:
        new S3RepositoryProvider(PROJECT, config).getName() == PROJECT

        where:
        PROJECT << [
            'my-bucket/simple',
            'my-bucket/org/team/project',
            'test-bucket/user/repo',
            'bucket/a/b/c/d'
        ]
    }
}
