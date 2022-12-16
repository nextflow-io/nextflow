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
    }

    def 'should set config' () {
        given:
        def OPTS = [
                debug:true,
                storageClass: 'FOO',
                storageKmsKeyId: 'key-1',
                storageEncryption: 'enc-2'
        ]

        when:
        def client = new AwsS3Config(OPTS)
        then:
        client.debug
        client.storageClass == 'FOO'
        client.storageKmsKeyId == 'key-1'
        client.storageEncryption == 'enc-2'
    }

    def 'should use legacy upload storage class' () {
        given:
        def OPTS = [
                uploadStorageClass: 'BAR',
        ]

        when:
        def client = new AwsS3Config(OPTS)
        then:
        client.storageClass == 'BAR'
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
}
