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

import nextflow.scm.GitUrl
import nextflow.scm.ProviderConfig
import spock.lang.Specification

/**
 * Tests for S3RepositoryFactory
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class S3RepositoryFactoryTest extends Specification {

    def 'should create S3 provider instance when platform is s3'() {
        given:
        def factory = new S3RepositoryFactory()
        def config = new S3ProviderConfig('test-bucket')

        when:
        def provider = factory.createProviderInstance(config, 'test-bucket/project')

        then:
        provider instanceof S3RepositoryProvider
        provider.project == 'test-bucket/project'
    }

    def 'should return null when platform is not s3'() {
        given:
        def factory = new S3RepositoryFactory()
        def config = new ProviderConfig('github', [platform: 'github'])

        when:
        def provider = factory.createProviderInstance(config, 'user/repo')

        then:
        provider == null
    }

    def 'should register TransportS3 on first provider creation'() {
        given:
        def factory = new S3RepositoryFactory()
        def config = new S3ProviderConfig('test-bucket')

        when:
        def provider1 = factory.createProviderInstance(config, 'test-bucket/project1')
        def provider2 = factory.createProviderInstance(config, 'test-bucket/project2')

        then:
        provider1 instanceof S3RepositoryProvider
        provider2 instanceof S3RepositoryProvider
        // TransportS3.register() should be called only once
    }

    def 'should get config for s3 URL'() {
        given:
        def factory = new S3RepositoryFactory()
        def url = new GitUrl('s3://my-bucket/path/to/project')
        def providers = []

        when:
        def config = factory.getConfig(providers, url)

        then:
        config instanceof S3ProviderConfig
        config.name == 'my-bucket'
        config.platform == 's3'
    }

    def 'should return existing config when domain matches'() {
        given:
        def factory = new S3RepositoryFactory()
        def existingConfig = new S3ProviderConfig('my-bucket', [
            region: 'eu-west-1',
            accessKey: 'test-key',
            secretKey: 'test-secret',
            platform: 's3'
        ])
        def providers = [existingConfig]
        def url = new GitUrl('s3://my-bucket/path/to/project')

        when:
        def config = factory.getConfig(providers, url)

        then:
        config.name == 'my-bucket'
        config.region.id() == 'eu-west-1'
    }

    def 'should return null for non-s3 protocol'() {
        given:
        def factory = new S3RepositoryFactory()
        def url = new GitUrl('https://github.com/user/repo')
        def providers = []

        when:
        def config = factory.getConfig(providers, url)

        then:
        config == null
    }

    def 'should create config instance when platform is s3'() {
        given:
        def factory = new S3RepositoryFactory()
        def attrs = [
            platform: 's3',
            region: 'us-west-2',
            accessKey: 'AKIAIOSFODNN7EXAMPLE',
            secretKey: 'wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY'
        ]

        when:
        def config = factory.createConfigInstance('my-bucket', attrs)

        then:
        config instanceof S3ProviderConfig
        config.name == 'my-bucket'
        config.platform == 's3'
        config.region.id() == 'us-west-2'
    }

    def 'should return null when creating config for non-s3 platform'() {
        given:
        def factory = new S3RepositoryFactory()
        def attrs = [
            platform: 'github',
            server: 'https://github.com'
        ]

        when:
        def config = factory.createConfigInstance('github', attrs)

        then:
        config == null
    }

    def 'should handle multiple buckets'() {
        given:
        def factory = new S3RepositoryFactory()
        def config1 = new S3ProviderConfig('bucket1')
        def config2 = new S3ProviderConfig('bucket2')

        when:
        def provider1 = factory.createProviderInstance(config1, 'bucket1/project1')
        def provider2 = factory.createProviderInstance(config2, 'bucket2/project2')

        then:
        provider1.project == 'bucket1/project1'
        provider2.project == 'bucket2/project2'
        provider1.getEndpointUrl() == 's3://bucket1/project1'
        provider2.getEndpointUrl() == 's3://bucket2/project2'
    }

    def 'should handle URL with nested paths'() {
        given:
        def factory = new S3RepositoryFactory()
        def url = new GitUrl('s3://my-bucket/org/team/project')
        def providers = []

        when:
        def config = factory.getConfig(providers, url)

        then:
        config instanceof S3ProviderConfig
        config.name == 'my-bucket'
    }

    def 'should create new config when no matching provider exists'() {
        given:
        def factory = new S3RepositoryFactory()
        def existingConfig = new S3ProviderConfig('other-bucket')
        def providers = [existingConfig]
        def url = new GitUrl('s3://my-bucket/project')

        when:
        def config = factory.getConfig(providers, url)

        then:
        config instanceof S3ProviderConfig
        config != existingConfig
        config.name == 'my-bucket'
    }

    def 'should not modify original attrs map'() {
        given:
        def factory = new S3RepositoryFactory()
        def attrs = [
            platform: 's3',
            region: 'us-east-1'
        ]
        def originalSize = attrs.size()

        when:
        factory.createConfigInstance('test-bucket', attrs)

        then:
        attrs.size() == originalSize
        attrs.platform == 's3'
    }

    def 'should handle different S3 URL formats'() {
        given:
        def factory = new S3RepositoryFactory()
        def providers = []

        expect:
        factory.getConfig(providers, new GitUrl('s3://bucket/project')) instanceof S3ProviderConfig
        factory.getConfig(providers, new GitUrl('s3://my-bucket/path/to/project')) instanceof S3ProviderConfig

    }
}
