/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.platform

import nextflow.SysEnv
import spock.lang.Specification

/**
 * Test PlatformHelper functionality
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
class PlatformHelperTest extends Specification {

    def 'should return Auth0 domain for known endpoints'() {
        expect:
        PlatformHelper.getAuthDomain(endpoint) == domain

        where:
        endpoint                              | domain
        'https://api.cloud.seqera.io'         | 'seqera.eu.auth0.com'
        'https://api.cloud.stage-seqera.io'   | 'seqera-stage.eu.auth0.com'
        'https://api.cloud.dev-seqera.io'     | 'seqera-development.eu.auth0.com'
        'https://enterprise.example.com'      | null
        'https://custom.endpoint.com'         | null
    }

    def 'should return Auth0 domain from environment variable'() {
        given:
        SysEnv.push(['TOWER_AUTH_DOMAIN': 'custom.auth0.com'])

        expect:
        PlatformHelper.getAuthDomain('https://custom.endpoint.com') == 'custom.auth0.com'

        cleanup:
        SysEnv.pop()
    }

    def 'should return Auth0 client ID for known endpoints'() {
        expect:
        PlatformHelper.getAuthClientId(endpoint) == clientId

        where:
        endpoint                              | clientId
        'https://api.cloud.seqera.io'         | 'FxCM8EJ76nNeHUDidSHkZfT8VtsrhHeL'
        'https://api.cloud.stage-seqera.io'   | '60cPDjI6YhoTPjyMTIBjGtxatSUwWswB'
        'https://api.cloud.dev-seqera.io'     | 'Ep2LhYiYmuV9hhz0dH6dbXVq0S7s7SWZ'
        'https://enterprise.example.com'      | null
        'https://custom.endpoint.com'         | null
    }

    def 'should return Auth0 client ID from environment variable'() {
        given:
        SysEnv.push(['TOWER_AUTH_CLIENT_ID': 'custom-client-id-123'])

        expect:
        PlatformHelper.getAuthClientId('https://custom.endpoint.com') == 'custom-client-id-123'

        cleanup:
        SysEnv.pop()
    }

    def 'should prioritize environment variable for Auth0 domain'() {
        given:
        SysEnv.push(['TOWER_AUTH_DOMAIN': 'override.auth0.com'])

        expect:
        PlatformHelper.getAuthDomain('https://api.cloud.seqera.io') == 'override.auth0.com'

        cleanup:
        SysEnv.pop()
    }

    def 'should prioritize environment variable for Auth0 client ID'() {
        given:
        SysEnv.push(['TOWER_AUTH_CLIENT_ID': 'override-client-id'])

        expect:
        PlatformHelper.getAuthClientId('https://api.cloud.seqera.io') == 'override-client-id'

        cleanup:
        SysEnv.pop()
    }
}
