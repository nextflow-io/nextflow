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

package nextflow.pixi

import java.nio.file.Paths

import nextflow.util.Duration
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Edmund Miller <edmund.miller@seqera.io>
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
        false       | [:]               | [NXF_PIXI_ENABLED: 'false']
        true        | [:]               | [NXF_PIXI_ENABLED: 'true']
        false       | [enabled: false]  | [NXF_PIXI_ENABLED: 'true']  // <-- config has priority
        true        | [enabled: true]   | [NXF_PIXI_ENABLED: 'true']
    }

    def 'should return create timeout'() {
        given:
        def CONFIG = [createTimeout: '30min']
        def pixi = new PixiConfig(CONFIG, [:])

        expect:
        pixi.createTimeout() == Duration.of('30min')
    }

    def 'should return null create timeout when not specified'() {
        given:
        def pixi = new PixiConfig([:], [:])

        expect:
        pixi.createTimeout() == null
    }

    def 'should return create options'() {
        given:
        def CONFIG = [createOptions: '--verbose --no-lock-update']
        def pixi = new PixiConfig(CONFIG, [:])

        expect:
        pixi.createOptions() == '--verbose --no-lock-update'
    }

    def 'should return null create options when not specified'() {
        given:
        def pixi = new PixiConfig([:], [:])

        expect:
        pixi.createOptions() == null
    }

    def 'should return cache directory'() {
        given:
        def CONFIG = [cacheDir: '/my/cache/dir']
        def pixi = new PixiConfig(CONFIG, [:])

        expect:
        pixi.cacheDir() == Paths.get('/my/cache/dir')
    }

    def 'should return null cache directory when not specified'() {
        given:
        def pixi = new PixiConfig([:], [:])

        expect:
        pixi.cacheDir() == null
    }

    def 'should handle boolean values for enabled flag'() {
        given:
        def pixi = new PixiConfig([enabled: true], [:])

        expect:
        pixi.isEnabled() == true

        when:
        pixi = new PixiConfig([enabled: false], [:])

        then:
        pixi.isEnabled() == false
    }

    def 'should handle string values for enabled flag'() {
        given:
        def pixi = new PixiConfig([enabled: 'true'], [:])

        expect:
        pixi.isEnabled() == true

        when:
        pixi = new PixiConfig([enabled: 'false'], [:])

        then:
        pixi.isEnabled() == false
    }

    def 'should inherit from LinkedHashMap'() {
        given:
        def CONFIG = [enabled: true, createTimeout: '10min', customOption: 'value']
        def pixi = new PixiConfig(CONFIG, [:])

        expect:
        pixi instanceof LinkedHashMap
        pixi.enabled == true
        pixi.createTimeout == '10min'
        pixi.customOption == 'value'
        pixi.size() == 3
    }
}
