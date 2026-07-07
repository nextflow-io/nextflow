/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.guix

import spock.lang.Specification
import spock.lang.Unroll

class GuixConfigTest extends Specification {

    @Unroll
    def 'should check enabled flag'() {
        given:
        def guix = new GuixConfig(CONFIG, ENV)
        expect:
        guix.isEnabled() == EXPECTED

        where:
        EXPECTED    | CONFIG            | ENV
        false       | [:]               | [:]
        false       | [enabled: false]  | [:]
        true        | [enabled: true]   | [:]
        and:
        false       | [:]               | [NXF_GUIX_ENABLED: 'false']
        true        | [:]               | [NXF_GUIX_ENABLED: 'true']
        false       | [enabled: false]  | [NXF_GUIX_ENABLED: 'true']  // <-- config has priority
        true        | [enabled: true]   | [NXF_GUIX_ENABLED: 'true']
    }

    def 'should check install options'() {
        given:
        def guix = new GuixConfig([installOptions: '--no-grafts'], [:])
        expect:
        guix.installOptions() == '--no-grafts'
    }

    def 'should have default create timeout'() {
        given:
        def guix = new GuixConfig([:], [:])
        expect:
        guix.createTimeout().toMillis() == 20 * 60 * 1000
    }
}
