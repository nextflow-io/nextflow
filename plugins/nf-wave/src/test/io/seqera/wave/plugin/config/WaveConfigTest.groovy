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

package io.seqera.wave.plugin.config

import io.seqera.wave.plugin.config.WaveConfig
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WaveConfigTest extends Specification {

    def 'should create empty opts' () {
        when:
        def opts = new WaveConfig([:])
        then:
        !opts.enabled()
        opts.endpoint() == 'http://localhost:9090'
    }

    def 'should create from env' () {
        given:
        def ENV = [WAVE_API_ENDPOINT: 'http://foo']
        when:
        def opts = new WaveConfig([:], ENV)
        then:
        !opts.enabled()
        opts.endpoint() == 'http://foo'
    }

    def 'should create config options' () {
        given:
        def ENV = [WAVE_API_ENDPOINT: 'http://foo']
        when:
        // config options have priority over sys env
        def opts = new WaveConfig([enabled:true, endpoint: 'http://localhost'], ENV)
        then:
        opts.enabled()
        opts.endpoint() == 'http://localhost'
    }

    def 'should remove ending slash' () {
        when:
        def opts = new WaveConfig([enabled:true, endpoint: 'http://localhost/v1//'])
        then:
        opts.enabled()
        opts.endpoint() == 'http://localhost/v1'
    }

    @Unroll
    def 'should add config urls'  () {
        when:
        def opts = new WaveConfig(OPTS, ENV)
        then:
        opts.containerConfigUrl() == EXPECTED

        where:
        OPTS                                                        | ENV                                                   | EXPECTED
        [:]                                                         | [:]                                                   | []
        [containerConfigUrl: 'http://foo.com']                      | [:]                                                   | [ new URL('http://foo.com')]
        [containerConfigUrl: 'http://foo.com']                      | [WAVE_CONTAINER_CONFIG_URL:'http://something.com']    | [ new URL('http://foo.com')]
        [:]                                                         | [WAVE_CONTAINER_CONFIG_URL:'http://something.com']    | [ new URL('http://something.com')]
        [containerConfigUrl: ['http://foo.com','https://bar.com']]  | [:]                                                   | [ new URL('http://foo.com'), new URL('https://bar.com')]
        [containerConfigUrl: ['http://foo.com','https://bar.com']]  | [WAVE_CONTAINER_CONFIG_URL:'http://boo.com']          | [ new URL('http://foo.com'), new URL('https://bar.com')]

    }

    def 'should get mamba config' () {
        when:
        def opts = new WaveConfig([:])
        then:
        opts.mambaOpts().from == 'mambaorg/micromamba:0.25.1'
        opts.mambaOpts().user == null

        when:
        opts = new WaveConfig([build:[mamba:[from:'mambaorg/foo:1', user:'hola']]])
        then:
        opts.mambaOpts().from == 'mambaorg/foo:1'
        opts.mambaOpts().user == 'hola'
        
    }

    def 'should set strategy' () {
        when:
        def opts = new WaveConfig([:])
        then:
        opts.strategy() == []

        when:
        opts = new WaveConfig([strategy:STRATEGY])
        then:
        opts.strategy() == EXPECTED

        where:
        STRATEGY                | EXPECTED
        null                    | []
        'dockerfile'            | ['dockerfile']
        'conda,container'       | ['conda','container']
        'conda , container'     | ['conda','container']
        ['conda','container']   | ['conda','container']
        [' conda',' container'] | ['conda','container']
    }

    def 'should fail to set strategy' () {
        when:
        new WaveConfig([strategy:['foo']])
        then:
        def e = thrown(IllegalArgumentException)
        e.message == "Invalid value for 'wave.strategy' configuration attribute - offending value: foo"
    }
}
