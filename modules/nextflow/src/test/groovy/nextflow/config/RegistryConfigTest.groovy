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

package nextflow.config

import spock.lang.Specification

/**
 * Tests for RegistryConfig
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class RegistryConfigTest extends Specification {

    def cleanup() {
        // Clean up any environment variables set during tests
        System.clearProperty('NXF_REGISTRY_TOKEN')
    }

    def 'should create config with default values'() {
        when:
        def config = new RegistryConfig()

        then:
        config.url == RegistryConfig.DEFAULT_REGISTRY_URL
        config.allUrls == [RegistryConfig.DEFAULT_REGISTRY_URL]
        config.getAuthToken(RegistryConfig.DEFAULT_REGISTRY_URL) == null
    }

    def 'should initialize with custom URL'() {
        given:
        def opts = [url: 'https://custom.registry.com']

        when:
        def config = new RegistryConfig(opts)

        then:
        config.url == 'https://custom.registry.com'
        config.allUrls == ['https://custom.registry.com']
    }

    def 'should initialize with multiple URLs'() {
        given:
        def opts = [
            urls: [
                'https://primary.registry.com',
                'https://fallback.registry.com'
            ]
        ]

        when:
        def config = new RegistryConfig(opts)

        then:
        config.allUrls == [
            'https://primary.registry.com',
            'https://fallback.registry.com'
        ]
    }

    def 'should prefer urls list over single url'() {
        given:
        def opts = [
            url: 'https://single.registry.com',
            urls: [
                'https://primary.registry.com',
                'https://fallback.registry.com'
            ]
        ]

        when:
        def config = new RegistryConfig(opts)

        then:
        config.allUrls == [
            'https://primary.registry.com',
            'https://fallback.registry.com'
        ]
    }

    def 'should fall back to default URL when no URL provided'() {
        given:
        def opts = [:]

        when:
        def config = new RegistryConfig(opts)

        then:
        config.url == RegistryConfig.DEFAULT_REGISTRY_URL
        config.allUrls == [RegistryConfig.DEFAULT_REGISTRY_URL]
    }

    def 'should initialize with authentication'() {
        given:
        def opts = [
            url: 'https://registry.com',
            auth: [
                'https://registry.com': 'token123'
            ]
        ]

        when:
        def config = new RegistryConfig(opts)

        then:
        config.getAuthToken('https://registry.com') == 'token123'
    }

    def 'should return null for unconfigured auth'() {
        given:
        def config = new RegistryConfig([url: 'https://registry.com'])

        when:
        def token = config.getAuthToken('https://registry.com')

        then:
        token == null
    }

    def 'should detect environment variable reference format'() {
        given:
        def opts = [
            url: 'https://registry.com',
            auth: [
                'https://registry.com': '${TEST_TOKEN_VAR}'
            ]
        ]
        def config = new RegistryConfig(opts)

        when:
        def token = config.getAuthToken('https://registry.com')

        then:
        token == '${TEST_TOKEN_VAR}'
        token.startsWith('${') && token.endsWith('}')
    }

    def 'should return literal token when not environment variable reference'() {
        given:
        def opts = [
            url: 'https://registry.com',
            auth: [
                'https://registry.com': 'literal-token'
            ]
        ]
        def config = new RegistryConfig(opts)

        when:
        def token = config.getAuthTokenResolved('https://registry.com')

        then:
        token == 'literal-token'
    }

    def 'should prefer configured auth over environment variable'() {
        given:
        def opts = [
            url: 'https://registry.com',
            auth: [
                'https://registry.com': 'config-token'
            ]
        ]
        def config = new RegistryConfig(opts)

        when:
        def token = config.getAuthTokenResolved('https://registry.com')

        then:
        // If NXF_REGISTRY_TOKEN env var exists, configured token takes precedence
        token == 'config-token'
    }

    def 'should check if auth is configured'() {
        given:
        def config = new RegistryConfig([
            url: 'https://registry.com',
            auth: [
                'https://registry.com': 'token123'
            ]
        ])

        expect:
        config.hasAuth('https://registry.com')
        !config.hasAuth('https://other-registry.com')
    }

    def 'should support multiple registry auths'() {
        given:
        def opts = [
            urls: [
                'https://primary.registry.com',
                'https://fallback.registry.com'
            ],
            auth: [
                'https://primary.registry.com': 'primary-token',
                'https://fallback.registry.com': 'fallback-token'
            ]
        ]
        def config = new RegistryConfig(opts)

        expect:
        config.getAuthToken('https://primary.registry.com') == 'primary-token'
        config.getAuthToken('https://fallback.registry.com') == 'fallback-token'
        config.hasAuth('https://primary.registry.com')
        config.hasAuth('https://fallback.registry.com')
    }

    def 'should handle null auth map'() {
        given:
        def config = new RegistryConfig([url: 'https://registry.com', auth: null])

        expect:
        config.getAuthToken('https://registry.com') == null
        !config.hasAuth('https://registry.com')
    }

    def 'should handle empty auth map'() {
        given:
        def config = new RegistryConfig([url: 'https://registry.com', auth: [:]])

        expect:
        config.getAuthToken('https://registry.com') == null
        !config.hasAuth('https://registry.com')
    }

    def 'should return null when environment variable is not set'() {
        given:
        def opts = [
            url: 'https://registry.com',
            auth: [
                'https://registry.com': '${NONEXISTENT_VAR}'
            ]
        ]
        def config = new RegistryConfig(opts)

        when:
        def token = config.getAuthTokenResolved('https://registry.com')

        then:
        token == null
    }

    def 'should handle complex URL patterns in auth keys'() {
        given:
        def opts = [
            auth: [
                'https://registry.com': 'token1',
                'https://registry.com:8080': 'token2',
                'http://localhost:3000': 'token3'
            ]
        ]
        def config = new RegistryConfig(opts)

        expect:
        config.getAuthToken('https://registry.com') == 'token1'
        config.getAuthToken('https://registry.com:8080') == 'token2'
        config.getAuthToken('http://localhost:3000') == 'token3'
    }

    def 'should handle empty urls list gracefully'() {
        given:
        def opts = [
            url: 'https://registry.com',
            urls: []
        ]

        when:
        def config = new RegistryConfig(opts)

        then:
        config.allUrls == ['https://registry.com']
    }

    def 'should use default URL when both url and urls are empty'() {
        given:
        def opts = [
            url: null,
            urls: []
        ]

        when:
        def config = new RegistryConfig(opts)

        then:
        config.allUrls == [RegistryConfig.DEFAULT_REGISTRY_URL]
    }

    def 'should correctly identify auth when configured'() {
        given:
        def config = new RegistryConfig([
            url: 'https://registry.com',
            auth: ['https://registry.com': 'token']
        ])

        expect:
        config.hasAuth('https://registry.com')
    }

    def 'should handle registry URL without trailing slash'() {
        given:
        def opts = [
            url: 'https://registry.com',
            auth: ['https://registry.com': 'token']
        ]
        def config = new RegistryConfig(opts)

        expect:
        config.getAuthToken('https://registry.com') == 'token'
    }

    def 'should handle registry URL with trailing slash'() {
        given:
        def opts = [
            url: 'https://registry.com/',
            auth: ['https://registry.com/': 'token']
        ]
        def config = new RegistryConfig(opts)

        expect:
        config.getAuthToken('https://registry.com/') == 'token'
    }

    def 'should preserve order of URLs in list'() {
        given:
        def opts = [
            urls: [
                'https://first.registry.com',
                'https://second.registry.com',
                'https://third.registry.com'
            ]
        ]

        when:
        def config = new RegistryConfig(opts)

        then:
        config.allUrls == [
            'https://first.registry.com',
            'https://second.registry.com',
            'https://third.registry.com'
        ]
    }
}
