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

import java.nio.file.Files

import nextflow.SysEnv
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsConfigTest extends Specification {

    @Unroll
    def 'should fetch aws creds'() {
        expect:
        AwsConfig.getAwsCredentials0(ENV, CREDS) == EXPECTED

        where:
        ENV                                                     | CREDS                                 | EXPECTED
        null                                                    | null                                  | null
        [AWS_ACCESS_KEY: 'x', AWS_SECRET_KEY: '222']            | null                                  | ['x','222']
        [AWS_ACCESS_KEY_ID: 'q', AWS_SECRET_ACCESS_KEY: '999']  | null                                  | ['q','999']
        [AWS_ACCESS_KEY: 'x', AWS_SECRET_KEY: '222',  AWS_ACCESS_KEY_ID: 'q', AWS_SECRET_ACCESS_KEY: '999']     | null | ['q','999']
        [AWS_ACCESS_KEY_ID: 'q', AWS_SECRET_ACCESS_KEY: '999']  | [accessKey: 'b', secretKey: '333']    | ['b','333']
        null                                                    | [accessKey: 'b', secretKey: '333']    | ['b','333']
        null                                                    | [accessKey: 'b']                      | null
        [AWS_ACCESS_KEY_ID: 'q', AWS_SECRET_ACCESS_KEY: '999']  | [accessKey: 'b', secretKey: '333']    | ['b','333']
    }

    def 'should lod creds from file'() {
        given:
        def file = Files.createTempFile('test','test')
        file.text = '''
            [default]
            aws_access_key_id = aaa
            aws_secret_access_key = bbbb
            '''

        expect:
        AwsConfig.getAwsCredentials0(null, null, [file]) == ['aaa','bbbb']
        and:
        AwsConfig.getAwsCredentials0([AWS_ACCESS_KEY: 'x', AWS_SECRET_KEY: '222'], null, [file]) == ['x','222']

        cleanup:
        file?.delete()
    }

    def 'should load creds with a profile'() {

        given:
        def file = Files.createTempFile('test','test')
        file.text = '''
            [default]
            aws_access_key_id = aaa
            aws_secret_access_key = bbbb
            
            [foo]
            aws_access_key_id = xxx
            aws_secret_access_key = yyy

            [bar]
            aws_access_key_id = ppp
            aws_secret_access_key = qqq
            aws_session_token = www
            '''

        expect:
        AwsConfig.getAwsCredentials0([AWS_PROFILE: 'foo'], null, [file]) == ['xxx','yyy']
        and:
        AwsConfig.getAwsCredentials0([AWS_DEFAULT_PROFILE: 'bar'], null, [file]) == ['ppp','qqq']

        cleanup:
        file?.delete()
    }

    def 'should get aws region'() {
        expect:
        AwsConfig.getAwsRegion([:], [:]) == null
        and:
        AwsConfig.getAwsRegion([:], [region:'eu-west-2']) == 'eu-west-2'
        and:
        // config has priority
        AwsConfig.getAwsRegion([AWS_DEFAULT_REGION: 'us-central-1'], [region:'eu-west-2']) == 'eu-west-2'
        and:
        AwsConfig.getAwsRegion([AWS_DEFAULT_REGION: 'us-central-1'], [:]) == 'us-central-1'
    }

    def 'should get aws region from aws file'() {
        given:
        def file = Files.createTempFile('test','test')
        file.text = '''
            [default]
            aws_access_key_id = aaa
            aws_secret_access_key = bbbb
            region = reg-something
            
            [foo]
            aws_access_key_id = xxx
            aws_secret_access_key = yyy
            region = reg-foo

            [bar]
            aws_access_key_id = xxx
            aws_secret_access_key = yyy
            aws_session_token = zzz
            '''

        expect:
        AwsConfig.getAwsRegion0([AWS_DEFAULT_REGION: 'us-central-1'], [:], file) == 'us-central-1'

        and:
        AwsConfig.getAwsRegion0([:], [:], file) == 'reg-something'

        and:
        AwsConfig.getAwsRegion0([:], [profile: 'foo'], file) == 'reg-foo'

        cleanup:
        file?.delete()
    }

    def 'should get aws credentials with file and profile from config' () {

        given:
        def file = Files.createTempFile('test','test')
        file.text = '''
            [default]
            aws_access_key_id = aaa
            aws_secret_access_key = bbbb
            
            [foo]
            aws_access_key_id = xxx
            aws_secret_access_key = yyy

            [bar]
            aws_access_key_id = ppp
            aws_secret_access_key = qqq
            aws_session_token = www
            '''

        expect:
        AwsConfig.getAwsCredentials0([:], [profile:'foo'], [file]) == ['xxx','yyy']
        AwsConfig.getAwsCredentials0([AWS_DEFAULT_PROFILE: 'bar'], [profile:'foo'], [file]) == ['xxx','yyy']
        AwsConfig.getAwsCredentials0([AWS_DEFAULT_PROFILE: 'bar'], [:], [file]) == ['ppp','qqq']
        AwsConfig.getAwsCredentials0([:], [:], [file]) == ['aaa','bbbb']

        cleanup:
        file?.delete()

    }

    def 'should get aws config' () {
        given:
        SysEnv.push(ENV)
        and:
        def config = new AwsConfig(CONFIG)

        expect:
        config.accessKey == ACCESS_KEY
        config.secretKey == SECRET_KEY
        config.profile == PROFILE
        config.region == REGION
        config.assumeRoleArn == ROLE

        cleanup:
        SysEnv.pop()

        where:
        ENV                 | CONFIG                                                                | ACCESS_KEY          | SECRET_KEY    | REGION            | PROFILE       | ROLE
        [:]                 | [accessKey: 'a', secretKey: 'b']                                      | 'a'                 | 'b'           | null              | 'default'     | null
        [:]                 | [accessKey: 'x', secretKey: 'y', region: 'eu-region-x']               | 'x'                 | 'y'           | 'eu-region-x'     | 'default'     | null
        [:]                 | [accessKey: 'p', secretKey: 'q', profile: 'hola']                     | 'p'                 | 'q'           | null              | 'hola'        | null
        [:]                 | [accessKey: 'a', secretKey: 'b', assumeRoleArn: 'role-x']             | 'a'                 | 'b'           | null              | 'default'     | 'role-x'
        and:
        [AWS_ACCESS_KEY_ID: 'k1', AWS_SECRET_ACCESS_KEY: 's1']  | [accessKey: 'a', secretKey: 'b']  | 'a'        | 'b'           | null              | 'default'     | null
        [AWS_ACCESS_KEY_ID: 'k1', AWS_SECRET_ACCESS_KEY: 's1']  | [:]                               | 'k1'              | 's1'          | null              | 'default'     | null
        [AWS_DEFAULT_REGION: 'eu-xyz']                          | [:]                               | null              | null          | 'eu-xyz'          | 'default'     | null
        [AWS_DEFAULT_PROFILE: 'my-profile']                     | [:]                               | null              | null          | null              | 'my-profile'  | null

    }

    @Unroll
    def 'should add max error retry' () {

        expect:
        AwsConfig.checkDefaultErrorRetry(SOURCE, ENV) == EXPECTED

        where:
        SOURCE                          | ENV                   | EXPECTED
        null                            | null                  | [max_error_retry: '5']
        [foo: 1]                        | [:]                   | [max_error_retry: '5', foo: 1]
        [foo: 1]                        | [AWS_MAX_ATTEMPTS:'3']| [max_error_retry: '3', foo: 1]
        [max_error_retry: '2', foo: 1]  | [:]                   | [max_error_retry: '2', foo: 1]
        [:]                             | [:]                   | [max_error_retry: '5']
    }
}
