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

import nextflow.config.spec.SpecNode
import nextflow.platform.AutoLabels
import nextflow.util.Duration
import spock.lang.Specification
import spock.lang.Unroll

/**
 * Unit tests for TowerConfig
 */
class TowerConfigTest extends Specification {

    def 'should use default endpoint when not specified'() {
        when:
        def config = new TowerConfig([:], [TOWER_API_ENDPOINT: 'https://example.com'])
        then:
        config.endpoint == 'https://example.com'

        when:
        config = new TowerConfig([:], [:])
        then:
        config.endpoint == 'https://api.cloud.seqera.io'
    }

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

    def 'should bound retries and timeouts for a best-effort lookup'() {
        given: 'options carrying the telemetry retry policy and long timeouts'
        def opts = [
            accessToken       : 'xyz',
            endpoint          : 'http://foo.com',
            httpConnectTimeout: Duration.of('60s'),
            httpReadTimeout   : Duration.of('60s') ]

        when:
        def config = TowerConfig.forLookup(opts, [:])

        then: 'the lookup gets a single attempt and short timeouts'
        config.retryPolicy.maxAttempts == 1
        config.httpConnectTimeout == TowerConfig.LOOKUP_TIMEOUT
        config.httpReadTimeout == TowerConfig.LOOKUP_TIMEOUT

        and: 'everything else is carried over unchanged'
        config.accessToken == 'xyz'
        config.endpoint == 'http://foo.com'

        and: 'the caller options are not modified'
        opts.httpConnectTimeout == Duration.of('60s')
        opts.httpReadTimeout == Duration.of('60s')
        !opts.containsKey('retryPolicy')
    }

    def 'should bound a lookup even when the config asks for many retries'() {
        when: 'the user configured an aggressive retry policy for telemetry'
        def config = TowerConfig.forLookup([accessToken: 'xyz', retryPolicy: [maxAttempts: 10]], [:])

        then: 'the lookup is still bounded to one attempt'
        config.retryPolicy.maxAttempts == 1
    }

    @Unroll
    def 'should keep the auto labels raw value for form: #VALUE'() {
        when:
        def config = new TowerConfig(OPTS, [:])

        then:
        // the value is carried over verbatim, because the runtime parses it with
        // `AutoLabels.parse` and needs to tell the accepted forms apart
        config.autoLabels == VALUE
        and:
        AutoLabels.parse(config.autoLabels, 'tower.autoLabels') == EXPECTED

        where:
        OPTS                                    | VALUE                        | EXPECTED
        [:]                                     | null                         | [] as Set
        [autoLabels: true]                      | true                         | AutoLabels.VALID_NAMES
        [autoLabels: false]                     | false                        | [] as Set
        [autoLabels: ['runName','projectName']] | ['runName','projectName']    | ['runName','projectName'] as Set
        [autoLabels: 'runName,projectName']     | 'runName,projectName'        | ['runName','projectName'] as Set
    }

    def 'should expose the auto labels option in the config spec'() {
        given:
        def root = SpecNode.Scope.of(TowerConfig, '')

        expect:
        // the option must be part of the spec, otherwise `ConfigValidator` flags it
        // as unrecognized and the config reference docs skip it
        root.getOption(['autoLabels']) != null
        and:
        // all the accepted forms are declared
        root.getOption(['autoLabels']).types().containsAll([Boolean, List, String])
        and:
        root.getOption(['autoLabels']).description().contains('`true`: include all available metadata labels')
        root.getOption(['autoLabels']).description().contains('`false` (default): disable')
        root.getOption(['autoLabels']).description().contains('a list or comma-separated string of short names')
        and:
        // every valid short name is documented
        AutoLabels.VALID_NAMES.every { root.getOption(['autoLabels']).description().contains("`${it}`") }
    }
}
