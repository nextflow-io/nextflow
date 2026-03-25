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

package io.seqera.tower.plugin

import nextflow.util.Duration
import spock.lang.Specification

/**
 * Unit tests for TowerConfig
 */
class TowerConfigTest extends Specification {

    def 'should use default timeout values when not specified'() {
        when:
        def config = new TowerConfig([:], [:])

        then:
        config.httpConnectTimeout == Duration.of('60s')
        config.httpReadTimeout == Duration.of('60s')
    }

    def 'should use provided connect timeout when specified'() {
        when:
        def config = new TowerConfig([httpConnectTimeout: Duration.of('30s')], [:])

        then:
        config.httpConnectTimeout == Duration.of('30s')
        config.httpReadTimeout == Duration.of('60s')
    }

    def 'should use provided read timeout when specified'() {
        when:
        def config = new TowerConfig([httpReadTimeout: Duration.of('120s')], [:])

        then:
        config.httpConnectTimeout == Duration.of('60s')
        config.httpReadTimeout == Duration.of('120s')
    }

    def 'should parse timeout from string value'() {
        when:
        def config = new TowerConfig([httpConnectTimeout: '5s', httpReadTimeout: '2m'], [:])

        then:
        config.httpConnectTimeout == Duration.of('5s')
        config.httpReadTimeout == Duration.of('2m')
    }
}
