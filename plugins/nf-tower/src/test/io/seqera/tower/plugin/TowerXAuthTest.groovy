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

import nextflow.util.ProxyConfig
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TowerXAuthTest extends Specification {

    def cleanup() {
        ProxyConfig.reset()
    }

    def 'should build a plain HttpClient without proxy when none is configured'() {
        when:
        def auth = new TowerXAuth('http://tower.example.com', 'tok', 'refresh')

        then:
        !auth.@httpClient.proxy().isPresent()
        !auth.@httpClient.authenticator().isPresent()
    }

    def 'should route the token-refresh HttpClient through the configured proxy'() {
        given:
        ProxyConfig.register(new ProxyConfig(protocol: 'https', host: 'proxy.example.com', port: '8080', username: 'foo', password: 'bar'))

        when:
        def auth = new TowerXAuth('http://tower.example.com', 'tok', 'refresh')

        then:
        auth.@httpClient.proxy().isPresent()
        auth.@httpClient.authenticator().isPresent()
    }

    def 'should route through the proxy without an authenticator when the proxy has no credentials'() {
        given:
        ProxyConfig.register(new ProxyConfig(protocol: 'https', host: 'proxy.example.com', port: '8080'))

        when:
        def auth = new TowerXAuth('http://tower.example.com', 'tok', 'refresh')

        then:
        auth.@httpClient.proxy().isPresent()
        !auth.@httpClient.authenticator().isPresent()
    }
}
