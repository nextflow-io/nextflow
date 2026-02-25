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

import nextflow.SysEnv
import spock.lang.Specification

/**
 * Tests for RegistryConfig
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class RegistryConfigTest extends Specification {

    def 'should create config with default values'() {
        when:
        def config = new RegistryConfig()

        then:
        config.url == RegistryConfig.DEFAULT_REGISTRY_URL
        config.allUrls == [RegistryConfig.DEFAULT_REGISTRY_URL]
        config.apiKey == null
    }

    def 'should initialize with custom URL as string'() {
        when:
        def config = new RegistryConfig([url: 'https://custom.registry.com'])

        then:
        config.url == 'https://custom.registry.com'
        config.allUrls == ['https://custom.registry.com']
    }

    def 'should initialize with URL as list'() {
        when:
        def config = new RegistryConfig([url: ['https://primary.registry.com', 'https://fallback.registry.com']])

        then:
        config.url == 'https://primary.registry.com'
        config.allUrls == ['https://primary.registry.com', 'https://fallback.registry.com']
    }

    def 'should use default URL when none provided'() {
        when:
        def config = new RegistryConfig([:])

        then:
        config.url == RegistryConfig.DEFAULT_REGISTRY_URL
        config.allUrls == [RegistryConfig.DEFAULT_REGISTRY_URL]
        config.apiKey == null
    }

    def 'should initialize with apiKey'() {
        when:
        def config = new RegistryConfig([url: 'https://registry.com', apiKey: 'token123'])

        then:
        config.apiKey == 'token123'
    }



    def 'should fall back to NXF_REGISTRY_TOKEN env var when apiKey not set'() {
        given:
        SysEnv.push([NXF_REGISTRY_TOKEN: 'env_var_token'])
        def config = new RegistryConfig([url: 'https://registry.com'])

        expect:
        // Without env var set, returns null
        config.apiKey == 'env_var_token'

        cleanup:
        SysEnv.pop()
    }

    def 'should preserve order of URLs in list'() {
        when:
        def config = new RegistryConfig([url: ['https://first.com', 'https://second.com', 'https://third.com']])

        then:
        config.allUrls == ['https://first.com', 'https://second.com', 'https://third.com']
        config.url == 'https://first.com'
    }

    def 'should handle empty list gracefully'() {
        when:
        def config = new RegistryConfig([url: []])

        then:
        config.allUrls == [RegistryConfig.DEFAULT_REGISTRY_URL]
        config.url == RegistryConfig.DEFAULT_REGISTRY_URL
    }
}