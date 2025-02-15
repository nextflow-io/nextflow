/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.cloud.aws.codecommit

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsCodeCommitProviderConfigTest extends Specification {

    def 'should create config' () {
        given:
        def HOST = 'git-codecommit.eu-west-1.amazonaws.com'
        when:
        def config = new AwsCodeCommitProviderConfig(HOST)
        then:
        config.name == 'codecommit'
        config.platform == 'codecommit'
        config.region == 'eu-west-1'
        config.domain == 'git-codecommit.eu-west-1.amazonaws.com'
        config.server == 'https://git-codecommit.eu-west-1.amazonaws.com'
        config.endpoint == 'https://git-codecommit.eu-west-1.amazonaws.com'

        expect:
        config.resolveProjectName('https://git-codecommit.eu-west-1.amazonaws.com/v1/repos/my-repo') == 'codecommit-eu-west-1/my-repo'
    }

}
