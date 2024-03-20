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
 */

package nextflow.cli

import nextflow.cli.v2.HubOptionsV2
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class HubOptionsTest extends Specification {

    def testUserV1() {

        when:
        def cmd = [:] as HubOptions.V1
        cmd.hubUserCli = credential
        then:
        cmd.getHubUser() == user
        cmd.getHubPassword() == password

        where:
        credential      | user     | password
        null            | null     | null
        'paolo'         | 'paolo'  | null
        'paolo:secret'  | 'paolo'  | 'secret'

    }

    def testUserV2() {

        when:
        def cmd = [:] as HubOptionsV2
        cmd.hubUserCli = credential
        then:
        cmd.getHubUser() == user
        cmd.getHubPassword() == password

        where:
        credential      | user     | password
        null            | null     | null
        'paolo'         | 'paolo'  | null
        'paolo:secret'  | 'paolo'  | 'secret'

    }

}
