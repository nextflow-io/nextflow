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

package nextflow.fusion

import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FusionConfigTest extends Specification {

    @Unroll
    def 'should create empty opts' () {
        when:
        def opts = new FusionConfig(OPTS, ENV)
        then:
        opts.enabled() == EXPECTED

        where:
        OPTS            | ENV           | EXPECTED
        [:]             | [:]           | false
        [enabled:false] | [:]           | false
        [enabled:true]  | [:]           | true
    }

    @Unroll
    def 'should create container config url' () {
        when:
        def opts = new FusionConfig(OPTS, ENV)
        then:
        opts.containerConfigUrl() == (EXPECTED ? new URL(EXPECTED) : null)

        where:
        OPTS                                    | ENV           | EXPECTED
        [:]                                     | [:]           | null
        [containerConfigUrl:'http://foo.com']   | [:]           | 'http://foo.com'
        [containerConfigUrl:'https://bar.com']  | [:]           | 'https://bar.com'
        [containerConfigUrl:'file:///some/file']| [:]           | 'file:///some/file'
        [:]                                     | [FUSION_CONTAINER_CONFIG_URL:'http://bar.com']           | 'http://bar.com'
        [containerConfigUrl:'http://foo.com']   | [FUSION_CONTAINER_CONFIG_URL:'http://bar.com']           | 'http://foo.com'

    }

    @Unroll
    def 'should get export aws key' () {
        expect:
        new FusionConfig(OPTS).exportStorageCredentials() == EXPECTED

        where:
        OPTS                            | EXPECTED
        [:]                             | false
        [exportAwsAccessKeys: false]    | false
        [exportAwsAccessKeys: true]     | true
    }

    @Unroll
    def 'should configure log level and output' () {
        given:
        def opts = new FusionConfig(OPTS)
        expect:
        opts.logLevel() == LEVEL
        opts.logOutput() == OUTPUT

        where:
        OPTS                            | LEVEL     | OUTPUT
        [:]                             | null      | null
        [logLevel: 'trace']             | 'trace'   | null
        [logOutput: 'stdout']           | null      | 'stdout'
    }

    def 'should configure cache size' () {
        given:
        def opts = new FusionConfig(OPTS)
        expect:
        opts.cacheSize() == SIZE

        where:
        OPTS                            | SIZE
        [:]                             | null
        [cacheSize: 100]                | MemoryUnit.of(100)
        [cacheSize: '100']              | MemoryUnit.of(100)
        [cacheSize: '100.MB']           | MemoryUnit.of('100.MB')
    }


    @Unroll
    def 'should configure tags' () {
        given:
        def opts = new FusionConfig(OPTS)
        expect:
        opts.tagsEnabled() == ENABLED
        opts.tagsPattern() == PATTERN

        where:
        OPTS                    | ENABLED   | PATTERN
        [:]                     | true      | FusionConfig.DEFAULT_TAGS
        [tags:true]             | true      | FusionConfig.DEFAULT_TAGS
        [tags:false]            | false     | null
        [tags:'[*.txt](x=1)']   | true      | '[*.txt](x=1)'

    }

    def 'should check privileged flag' () {
        given:
        def opts = new FusionConfig(OPTS)
        expect:
        opts.privileged() == EXPECTED

        where:
        OPTS                    | EXPECTED
        [:]                     | true
        [privileged:true]       | true
        [privileged:false]      | false
    }

    @Unroll
    def 'should parse fusion version' () {
        expect:
        new FusionConfig([:]).retrieveFusionVersion(FUSION_URL) == EXPECTED
        where:
        FUSION_URL                              | EXPECTED
        FusionConfig.DEFAULT_FUSION_AMD64_URL   | '2.4'
        FusionConfig.DEFAULT_FUSION_ARM64_URL   | '2.4'
        'https://foo.com/releases/v3.0-amd.json'| '3.0'
    }

    def 'should get version version from config' () {
        expect:
        new FusionConfig([containerConfigUrl:FUSION_URL, enabled:ENABLED]).version() == EXPECTED
        where:
        FUSION_URL                                      | ENABLED  | EXPECTED
        null                                            | false    | null
        null                                            | true     | '2.4'
        'https://foo.com/releases/v4.0-amd64.json'      | true     | '4.0'
        'https://foo.com/releases/v4.0.1-amd64.json'    | true     | '4.0.1'
    }
}
