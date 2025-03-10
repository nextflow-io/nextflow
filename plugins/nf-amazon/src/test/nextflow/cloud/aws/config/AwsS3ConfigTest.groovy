/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.cloud.aws.config

import com.amazonaws.services.s3.model.CannedAccessControlList
import nextflow.SysEnv
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsS3ConfigTest extends Specification {

    def 'should get default config' () {
        when:
        def client = new AwsS3Config([:])
        then:
        !client.storageClass
        !client.storageKmsKeyId
        !client.storageEncryption
        !client.debug
        !client.s3Acl
        !client.pathStyleAccess
        !client.anonymous
        !client.isCustomEndpoint()
    }

    def 'should set config' () {
        given:
        def OPTS = [
                debug:true,
                storageClass: 'STANDARD',
                storageKmsKeyId: 'key-1',
                storageEncryption: 'AES256',
                s3Acl: 'public-read',
                s3PathStyleAccess: true,
                anonymous: true
        ]

        when:
        def client = new AwsS3Config(OPTS)
        then:
        client.debug
        client.storageClass == 'STANDARD'
        client.storageKmsKeyId == 'key-1'
        client.storageEncryption == 'AES256'
        client.s3Acl == CannedAccessControlList.PublicRead
        client.pathStyleAccess
        client.anonymous
    }

    def 'should use legacy upload storage class' () {
        given:
        def OPTS = [
                uploadStorageClass: 'STANDARD_IA',
        ]

        when:
        def client = new AwsS3Config(OPTS)
        then:
        client.storageClass == 'STANDARD_IA'
    }

    @Unroll
    def 'should get aws s3 endpoint' () {
        given:
        SysEnv.push(ENV)

        when:
        def config = new AwsS3Config(CONFIG)
        then:
        config.endpoint == EXPECTED

        cleanup:
        SysEnv.pop()

        where:
        ENV                             | CONFIG                        | EXPECTED
        [:]                             | [:]                           | null
        [AWS_S3_ENDPOINT: 'http://foo'] | [:]                           | 'http://foo'
        [:]                             | [endpoint: 'http://bar']      | 'http://bar'
        [AWS_S3_ENDPOINT: 'http://foo'] | [endpoint: 'http://bar']      | 'http://bar'  // <-- config should have priority
    }

    @Unroll
    def 'should fail with invalid endpoint protocol' () {
        when:
        new AwsS3Config(CONFIG)
        then:
        def e = thrown(IllegalArgumentException)
        e.message == EXPECTED

        where:
        CONFIG                          | EXPECTED
        [endpoint: 'bar.com']           |  "S3 endpoint must begin with http:// or https:// prefix - offending value: 'bar.com'"
        [endpoint: 'ftp://bar.com']     |  "S3 endpoint must begin with http:// or https:// prefix - offending value: 'ftp://bar.com'"

    }

    def 'should get s3 legacy properties' () {
        given:
        SysEnv.push([:])

        when:
        def config = new AwsConfig([client:[uploadMaxThreads: 5, uploadChunkSize: 1000, uploadStorageClass: 'STANDARD']])
        def env = config.getS3LegacyProperties()
        then:
        env.upload_storage_class == 'STANDARD'
        env.upload_chunk_size == '1000'
        env.upload_max_threads == '5'
        env.max_error_retry == '5'  // <-- default to 5

        when:
        config = new AwsConfig([client:[uploadMaxThreads: 10, maxErrorRetry: 20, uploadStorageClass: 'ONEZONE_IA']])
        env = config.getS3LegacyProperties()

        then:
        env.upload_storage_class == 'ONEZONE_IA'
        env.upload_max_threads == '10'
        env.max_error_retry == '20'

        cleanup:
        SysEnv.pop()

    }

    @Unroll
    def 'should check is custom endpoint' () {
        given:
        def config = new AwsS3Config(CONFIG)

        expect:
        config.isCustomEndpoint() == EXPECTED

        where:
        EXPECTED    | CONFIG
        false       | [:]
        false       | [endpoint: 'https://s3.us-east-2.amazonaws.com']
        true        | [endpoint: 'https://foo.com']
        // consider AWS china as custom ednpoint
        // see https://github.com/nextflow-io/nextflow/issues/5836
        true        | [endpoint: 'https://xxxx.s3.cn-north-1.vpce.amazonaws.com.cn']
    }
}
