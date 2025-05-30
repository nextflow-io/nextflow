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

package nextflow.pixi

import spock.lang.Specification
import spock.lang.Unroll

/**
 * Tests for PixiConfig class
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PixiConfigTest extends Specification {

    @Unroll
    def 'should check enabled flag'() {
        given:
        def pixi = new PixiConfig(CONFIG, ENV)
        expect:
        pixi.isEnabled() == EXPECTED

        where:
        EXPECTED    | CONFIG            | ENV
        false       | [:]               | [:]
        false       | [enabled: false]  | [:]
        true        | [enabled: true]   | [:]
        and:
        false       | [:]               | [NXF_PIXI_ENABLED: false]
        true        | [:]               | [NXF_PIXI_ENABLED: true]
        false       | [enabled: false]  | [NXF_PIXI_ENABLED: true]  // <-- config has priority
        true        | [enabled: true]   | [NXF_PIXI_ENABLED: true]
    }

    def 'should create config with defaults'() {
        when:
        def pixi = new PixiConfig([:])
        then:
        !pixi.isEnabled()
    }

    def 'should create config with enabled true'() {
        when:
        def pixi = new PixiConfig([enabled: true])
        then:
        pixi.isEnabled()
    }

    def 'should create config with environment variable'() {
        when:
        def pixi = new PixiConfig([:], [NXF_PIXI_ENABLED: 'true'])
        then:
        pixi.isEnabled()
    }

    def 'should prioritize config over environment'() {
        when:
        def pixi = new PixiConfig([enabled: false], [NXF_PIXI_ENABLED: 'true'])
        then:
        !pixi.isEnabled()
    }
}
