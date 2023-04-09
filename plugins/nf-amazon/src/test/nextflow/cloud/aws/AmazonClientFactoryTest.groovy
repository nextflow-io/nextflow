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


import nextflow.cloud.aws.config.AwsConfig
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AmazonClientFactoryTest extends Specification {

    def 'should create factory' () {
        when:
        def factory = new AmazonClientFactory(new AwsConfig(accessKey: 'foo', secretKey: 'bar', region:'xyz', profile:'my-profile'))
        then:
        factory.accessKey() == 'foo'
        factory.secretKey() == 'bar'
        factory.region() == 'xyz'
        factory.profile() == 'my-profile'
    }

}
