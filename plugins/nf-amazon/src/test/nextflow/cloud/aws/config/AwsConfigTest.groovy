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
        config.credentials == (ACCESS_KEY && SECRET_KEY ? [ACCESS_KEY, SECRET_KEY] : [])

        cleanup:
        SysEnv.pop()

        where:
        ENV                 | CONFIG                                                                | ACCESS_KEY          | SECRET_KEY    | REGION            | PROFILE
        [:]                 | [accessKey: 'a', secretKey: 'b']                                      | 'a'                 | 'b'           | null              | null
        [:]                 | [accessKey: 'x', secretKey: 'y', region: 'eu-region-x']               | 'x'                 | 'y'           | 'eu-region-x'     | null
        [:]                 | [accessKey: 'p', secretKey: 'q', profile: 'hola']                     | 'p'                 | 'q'           | null              | 'hola'
        and:
        [AWS_DEFAULT_REGION: 'eu-xyz']                          | [:]                               | null              | null          | 'eu-xyz'          | null
        [AWS_DEFAULT_PROFILE: 'my-profile']                     | [:]                               | null              | null          | null              | 'my-profile'

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
