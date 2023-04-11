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

package nextflow.cloud.aws

import nextflow.SysEnv
import nextflow.cloud.aws.config.AwsConfig
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsClientFactoryTest extends Specification {

    def 'should create factory' () {
        given:
        SysEnv.push([:])
        when:
        def factory = new AwsClientFactory(new AwsConfig(accessKey: 'foo', secretKey: 'bar', region:'xyz', profile:'my-profile'))
        then:
        factory.accessKey() == 'foo'
        factory.secretKey() == 'bar'
        factory.region() == 'xyz'
        factory.profile() == 'my-profile'

        cleanup:
        SysEnv.pop()
    }

    def 'should create factory using environment' () {
        given:
        SysEnv.push([AWS_REGION:'eu-foo-1', AWS_PROFILE: 'profile-x'])
        when:
        def factory = new AwsClientFactory(new AwsConfig([:]))
        then:
        factory.accessKey() == null
        factory.secretKey() == null
        factory.region() == 'eu-foo-1'
        factory.profile() == 'profile-x'

        cleanup:
        SysEnv.pop()
    }
}
