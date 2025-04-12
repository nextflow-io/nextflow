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

package nextflow.script

import nextflow.Session
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FusionMetaTest extends Specification {

    @Unroll
    def 'should get fusion meta' () {
        given:
        def session = Mock(Session) { getConfig()>>OPTS }

        expect:
        new FusionMetadata(session).enabled == EXPECTED_ENABLED
        new FusionMetadata(session).version == EXPECTED_VERSION

        where:
        OPTS                        | EXPECTED_ENABLED  | EXPECTED_VERSION
        [:]                         | false             | null
        [fusion:[enabled:false]]    | false             | null
        [fusion:[enabled:true]]     | true              | '2.4'
        [fusion:[enabled:true, containerConfigUrl: 'https://foo.io/releases/v3.0-amd64.json']]     | true    | '3.0'
    }

}
