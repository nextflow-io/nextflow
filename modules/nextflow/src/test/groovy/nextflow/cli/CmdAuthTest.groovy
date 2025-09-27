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

package nextflow.cli

import nextflow.exception.AbortOperationException
import org.junit.Rule
import spock.lang.Specification
import spock.lang.TempDir
import test.OutputCapture

import java.nio.file.Files
import java.nio.file.Path

/**
 * Test CmdAuth functionality
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
class CmdAuthTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @TempDir
    Path tempDir

    def 'should have correct name'() {
        given:
        def cmd = new CmdAuth()

        expect:
        cmd.getName() == 'auth'
    }

    def 'should define correct constants'() {
        expect:
        CmdAuth.SEQERA_ENDPOINTS.prod == 'https://api.cloud.seqera.io'
        CmdAuth.SEQERA_ENDPOINTS.stage == 'https://api.cloud.stage-seqera.io'
        CmdAuth.SEQERA_ENDPOINTS.dev == 'https://api.cloud.dev-seqera.io'
        CmdAuth.API_TIMEOUT_MS == 10000
        CmdAuth.AUTH_POLL_TIMEOUT_RETRIES == 60
        CmdAuth.AUTH_POLL_INTERVAL_SECONDS == 5
        CmdAuth.WORKSPACE_SELECTION_THRESHOLD == 8
    }

    def 'should show usage when no args provided'() {
        given:
        def cmd = new CmdAuth()

        when:
        cmd.run()

        then:
        def output = capture.toString()
        output.contains('Manage Seqera Platform authentication')
        output.contains('Usage: nextflow auth <sub-command> [options]')
        output.contains('Commands:')
        output.contains('login')
        output.contains('logout')
        output.contains('config')
        output.contains('status')
    }

    def 'should show specific command usage'() {
        given:
        def cmd = new CmdAuth()
        cmd.args = ['login']

        when:
        cmd.usage()

        then:
        def output = capture.toString()
        output.contains('Authenticate with Seqera Platform')
        output.contains('Usage: nextflow auth login')
        output.contains('-u, -url <endpoint>')
    }

    def 'should throw error for unknown command'() {
        given:
        def cmd = new CmdAuth()
        cmd.args = ['unknown']

        when:
        cmd.run()

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('Unknown auth sub-command: unknown')
    }

    def 'should suggest closest command for typos'() {
        given:
        def cmd = new CmdAuth()

        when:
        cmd.getCmd(['loginn'])

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('Unknown auth sub-command: loginn')
        ex.message.contains('Did you mean one of these?')
        ex.message.contains('login')
    }

    def 'should identify cloud endpoints correctly'() {
        given:
        def cmd = new CmdAuth()

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
        def cmd = new CmdAuth()

        expect:
        cmd.isCloudEndpoint('https://api.cloud.seqera.io') == true
        cmd.isCloudEndpoint('https://api.cloud.stage-seqera.io') == true
        cmd.isCloudEndpoint('https://enterprise.example.com') == false
    }

    def 'should validate argument count correctly'() {
        given:
        def cmd = new CmdAuth()

        when:
        cmd.validateArgumentCount(['extra'], 'test')

        then:
        def ex = thrown(AbortOperationException)
        ex.message == 'Too many arguments for test command'

        when:
        cmd.validateArgumentCount([], 'test')

        then:
        noExceptionThrown()
    }

    def 'should read config correctly'() {
        given:
        def cmd = new CmdAuth()

        expect:
        // readConfig method should return a Map
        cmd.readConfig() instanceof Map
    }

    def 'should clean tower config from existing content'() {
        given:
        def cmd = new CmdAuth()
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
        def cmd = new CmdAuth()
        def config = [
            'tower.accessToken': 'test-token',
            'tower.enabled': true
        ]

        when:
        cmd.writeConfig(config, null)

        then:
        noExceptionThrown()
    }

    def 'login command should validate too many arguments'() {
        given:
        def cmd = new CmdAuth()
        cmd.args = ['login', 'extra']

        when:
        cmd.run()

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('Too many arguments for login command')
    }

    def 'logout command should validate too many arguments'() {
        given:
        def cmd = new CmdAuth()
        cmd.args = ['logout', 'extra']

        when:
        cmd.run()

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('Too many arguments for logout command')
    }

    def 'config command should validate too many arguments'() {
        given:
        def cmd = new CmdAuth()
        cmd.args = ['config', 'extra']

        when:
        cmd.run()

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('Too many arguments for config command')
    }

    def 'status command should validate too many arguments'() {
        given:
        def cmd = new CmdAuth()
        cmd.args = ['status', 'extra']

        when:
        cmd.run()

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('Too many arguments for status command')
    }

    def 'login command should use provided API URL'() {
        given:
        def cmd = new CmdAuth()
        cmd.args = ['login']
        cmd.apiUrl = 'https://api.example.com'

        when:
        def loginCmd = cmd.getCmd(cmd.args)

        then:
        loginCmd instanceof CmdAuth.LoginCmd

        when:
        cmd.run()

        then:
        loginCmd.apiUrl == 'https://api.example.com'
    }

    def 'should have all required subcommands'() {
        given:
        def cmd = new CmdAuth()

        expect:
        cmd.commands.size() == 4
        cmd.commands.find { it.name == 'login' } != null
        cmd.commands.find { it.name == 'logout' } != null
        cmd.commands.find { it.name == 'config' } != null
        cmd.commands.find { it.name == 'status' } != null
    }

    def 'should create HTTP connection with correct properties'() {
        given:
        def cmd = new CmdAuth()

        when:
        def connection = cmd.createHttpConnection('https://example.com', 'GET', 'test-token')

        then:
        connection.requestMethod == 'GET'
        connection.connectTimeout == CmdAuth.API_TIMEOUT_MS
        connection.readTimeout == CmdAuth.API_TIMEOUT_MS

        when:
        def connectionNoAuth = cmd.createHttpConnection('https://example.com', 'POST')

        then:
        connectionNoAuth.requestMethod == 'POST'
        connectionNoAuth.connectTimeout == CmdAuth.API_TIMEOUT_MS
        connectionNoAuth.readTimeout == CmdAuth.API_TIMEOUT_MS
    }
}