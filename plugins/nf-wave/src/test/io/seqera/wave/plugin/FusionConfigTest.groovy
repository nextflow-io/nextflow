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

package io.seqera.wave.plugin

import io.seqera.wave.plugin.config.FusionConfig
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
        [:]                                     | [FUSION_CONTAINER_CONFIG_URL:'http://bar.com']           | 'http://bar.com'
        [containerConfigUrl:'http://foo.com']   | [FUSION_CONTAINER_CONFIG_URL:'http://bar.com']           | 'http://foo.com'

    }


}
