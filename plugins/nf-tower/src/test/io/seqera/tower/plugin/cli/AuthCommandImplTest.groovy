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

package io.seqera.tower.plugin.cli

import nextflow.exception.AbortOperationException
import org.junit.Rule
import spock.lang.Specification
import spock.lang.TempDir
import test.OutputCapture

import java.nio.file.Path

/**
 * Test CmdAuth functionality
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
class AuthCommandImplTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @TempDir
    Path tempDir





    def 'should identify cloud endpoints correctly'() {
        given:
        def cmd = new AuthCommandImpl()

        expect:
        cmd.getCloudEndpointInfo('https://api.cloud.seqera.io').isCloud == true
        cmd.getCloudEndpointInfo('https://api.cloud.seqera.io').environment == 'prod'
        cmd.getCloudEndpointInfo('https://api.cloud.stage-seqera.io').isCloud == true
        cmd.getCloudEndpointInfo('https://api.cloud.stage-seqera.io').environment == 'stage'
        cmd.getCloudEndpointInfo('https://api.cloud.dev-seqera.io').isCloud == true
        cmd.getCloudEndpointInfo('https://api.cloud.dev-seqera.io').environment == 'dev'
        cmd.getCloudEndpointInfo('https://cloud.seqera.io/api').isCloud == true
        cmd.getCloudEndpointInfo('https://cloud.seqera.io/api').environment == 'prod'
        cmd.getCloudEndpointInfo('https://enterprise.example.com').isCloud == false
        cmd.getCloudEndpointInfo('https://enterprise.example.com').environment == null
    }

    def 'should identify cloud endpoint from URL'() {
        given:
        def cmd = new AuthCommandImpl()

        expect:
        cmd.isCloudEndpoint('https://api.cloud.seqera.io') == true
        cmd.isCloudEndpoint('https://api.cloud.stage-seqera.io') == true
        cmd.isCloudEndpoint('https://enterprise.example.com') == false
    }

    def 'should read config correctly'() {
        given:
        def cmd = new AuthCommandImpl()

        expect:
        // readConfig method should return a Map
        cmd.readConfig() instanceof Map
    }

    def 'should clean tower config from existing content'() {
        given:
        def cmd = new AuthCommandImpl()
        def content = '''
// Some other config
process {
    executor = 'local'
}

// Seqera Platform configuration
tower {
    accessToken = 'old-token'
    enabled = true
}

tower.endpoint = 'old-endpoint'

// More config
params.test = true
'''

        when:
        def cleaned = cmd.cleanTowerConfig(content)

        then:
        !cleaned.contains('tower {')
        !cleaned.contains('accessToken = \'old-token\'')
        !cleaned.contains('tower.endpoint')
        !cleaned.contains('Seqera Platform configuration')
        cleaned.contains('process {')
        cleaned.contains('params.test = true')
    }

    def 'should handle config writing'() {
        given:
        def cmd = new AuthCommandImpl()
        def config = [
            'tower.accessToken': 'test-token',
            'tower.enabled': true
        ]

        when:
        cmd.writeConfig(config, null)

        then:
        noExceptionThrown()
    }

    def 'should create HTTP connection with correct properties'() {
        given:
        def cmd = new AuthCommandImpl()

        when:
        def connection = cmd.createHttpConnection('https://example.com', 'GET', 'test-token')

        then:
        connection.requestMethod == 'GET'
        connection.connectTimeout == AuthCommandImpl.API_TIMEOUT_MS
        connection.readTimeout == AuthCommandImpl.API_TIMEOUT_MS

        when:
        def connectionNoAuth = cmd.createHttpConnection('https://example.com', 'POST')

        then:
        connectionNoAuth.requestMethod == 'POST'
        connectionNoAuth.connectTimeout == AuthCommandImpl.API_TIMEOUT_MS
        connectionNoAuth.readTimeout == AuthCommandImpl.API_TIMEOUT_MS
    }
}