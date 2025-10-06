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

import io.seqera.http.HxClient
import nextflow.Const
import nextflow.SysEnv
import org.junit.Rule
import spock.lang.Specification
import spock.lang.TempDir
import test.OutputCapture

import java.net.http.HttpResponse
import java.nio.file.Files
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
        cmd.getCloudEndpointInfo('https://api.cloud.seqera.io').auth.domain == 'seqera.eu.auth0.com'
        cmd.getCloudEndpointInfo('https://api.cloud.stage-seqera.io').isCloud == true
        cmd.getCloudEndpointInfo('https://api.cloud.stage-seqera.io').auth.domain == 'seqera-stage.eu.auth0.com'
        cmd.getCloudEndpointInfo('https://api.cloud.dev-seqera.io').isCloud == true
        cmd.getCloudEndpointInfo('https://api.cloud.dev-seqera.io').auth.domain == 'seqera-development.eu.auth0.com'
        cmd.getCloudEndpointInfo('https://cloud.seqera.io/api').isCloud == true
        cmd.getCloudEndpointInfo('https://cloud.seqera.io/api').auth.domain == 'seqera.eu.auth0.com'
        cmd.getCloudEndpointInfo('https://enterprise.example.com').isCloud == false
        cmd.getCloudEndpointInfo('https://enterprise.example.com').auth == null
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

    def 'should handle config writing'() {
        given:
        def cmd = new AuthCommandImpl()
        def config = [
            'tower.accessToken': 'test-token',
            'tower.enabled'    : true
        ]

        when:
        cmd.writeConfig(config, null)

        then:
        noExceptionThrown()
    }


    def 'should normalize API URL correctly'() {
        given:
        def cmd = new AuthCommandImpl()

        expect:
        cmd.normalizeApiUrl(null) == 'https://api.cloud.seqera.io'
        cmd.normalizeApiUrl('') == 'https://api.cloud.seqera.io'
        cmd.normalizeApiUrl('https://example.com') == 'https://example.com'
        cmd.normalizeApiUrl('http://example.com') == 'http://example.com'
        cmd.normalizeApiUrl('example.com') == 'https://example.com'
        cmd.normalizeApiUrl('api.cloud.seqera.io') == 'https://api.cloud.seqera.io'
    }

    def 'should decode token ID correctly'() {
        given:
        def cmd = new AuthCommandImpl()
        def tokenData = '{"tid":"token-123","exp":1234567890}'
        def encodedToken = Base64.encoder.encodeToString(tokenData.bytes)

        when:
        def tokenId = cmd.decodeTokenId(encodedToken)

        then:
        tokenId == 'token-123'
    }

    def 'should fail to decode token without ID'() {
        given:
        def cmd = new AuthCommandImpl()
        def tokenData = '{"exp":1234567890}'
        def encodedToken = Base64.encoder.encodeToString(tokenData.bytes)

        when:
        cmd.decodeTokenId(encodedToken)

        then:
        thrown(RuntimeException)
    }

    def 'should fail to decode invalid token'() {
        given:
        def cmd = new AuthCommandImpl()
        def invalidToken = 'not-a-valid-base64-token'

        when:
        cmd.decodeTokenId(invalidToken)

        then:
        thrown(RuntimeException)
    }

    def 'should remove includeConfig line correctly'() {
        given:
        def cmd = new AuthCommandImpl()
        def content = """
// Some config
param1 = 'value1'

includeConfig 'seqera_auth.config'

param2 = 'value2'
"""

        when:
        def result = cmd.removeIncludeConfigLine(content)

        then:
        !result.contains("includeConfig 'seqera_auth.config'")
        result.contains('param1 = \'value1\'')
        result.contains('param2 = \'value2\'')
    }

    def 'should handle content without includeConfig line'() {
        given:
        def cmd = new AuthCommandImpl()
        def content = """
// Some config
param1 = 'value1'
param2 = 'value2'"""

        when:
        def result = cmd.removeIncludeConfigLine(content)

        then:
        result == content
    }

    def 'should get config file path'() {
        given:
        def cmd = new AuthCommandImpl()

        when:
        def configFile = cmd.getConfigFile()

        then:
        configFile == Const.APP_HOME_DIR.resolve('config')
    }

    def 'should get auth file path'() {
        given:
        def cmd = new AuthCommandImpl()

        when:
        def authFile = cmd.getAuthFile()

        then:
        authFile == Const.APP_HOME_DIR.resolve('seqera_auth.config')
    }

    def 'should read empty auth file'() {
        given:
        def cmd = new AuthCommandImpl()

        when:
        def config = cmd.readAuthFile()

        then:
        config instanceof Map
        config.isEmpty() || config.size() >= 0
    }

    def 'should write config to seqera_auth.config file'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def authFile = tempDir.resolve('seqera_auth.config')
        def configFile = tempDir.resolve('config')

        cmd.getAuthFile() >> authFile
        cmd.getConfigFile() >> configFile

        def config = [
            'tower.accessToken': 'test-token-123',
            'tower.endpoint'   : 'https://api.cloud.seqera.io',
            'tower.enabled'    : true
        ]

        when:
        cmd.writeConfig(config, null)

        then:
        Files.exists(authFile)
        def authContent = Files.readString(authFile)
        authContent.contains('accessToken = \'test-token-123\'')
        authContent.contains('endpoint = \'https://api.cloud.seqera.io\'')
        authContent.contains('enabled = true')
    }

    def 'should write config with workspace metadata'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def authFile = tempDir.resolve('seqera_auth.config')
        def configFile = tempDir.resolve('config')

        cmd.getAuthFile() >> authFile
        cmd.getConfigFile() >> configFile

        def config = [
            'tower.accessToken': 'test-token-123',
            'tower.endpoint'   : 'https://api.cloud.seqera.io',
            'tower.workspaceId': '12345'
        ]
        def metadata = [
            orgName          : 'TestOrg',
            workspaceName    : 'TestWorkspace',
            workspaceFullName: 'test-org/test-workspace'
        ]

        when:
        cmd.writeConfig(config, metadata)

        then:
        Files.exists(authFile)
        def authContent = Files.readString(authFile)
        authContent.contains('workspaceId = \'12345\'')
        authContent.contains('// TestOrg / TestWorkspace [test-org/test-workspace]')
    }

    def 'should add includeConfig line to main config file'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def authFile = tempDir.resolve('seqera_auth.config')
        def configFile = tempDir.resolve('config')

        cmd.getAuthFile() >> authFile
        cmd.getConfigFile() >> configFile

        // Create initial config file
        Files.writeString(configFile, "// Existing config\nparam1 = 'value1'\n")

        def config = [
            'tower.accessToken': 'test-token-123'
        ]

        when:
        cmd.writeConfig(config, null)

        then:
        Files.exists(configFile)
        def configContent = Files.readString(configFile)
        configContent.contains("includeConfig 'seqera_auth.config'")
        configContent.contains('param1 = \'value1\'')
    }

    def 'should not duplicate includeConfig line'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def authFile = tempDir.resolve('seqera_auth.config')
        def configFile = tempDir.resolve('config')

        cmd.getAuthFile() >> authFile
        cmd.getConfigFile() >> configFile

        // Create config file with existing includeConfig line
        Files.writeString(configFile, "// Config\nincludeConfig 'seqera_auth.config'\nparam1 = 'value1'\n")

        def config = [
            'tower.accessToken': 'test-token-123'
        ]

        when:
        cmd.writeConfig(config, null)

        then:
        Files.exists(configFile)
        def configContent = Files.readString(configFile)
        configContent.count("includeConfig 'seqera_auth.config'") == 1
    }

    def 'should get current workspace name when workspace exists'() {
        given:
        def cmd = new AuthCommandImpl()
        def workspaces = [
            [workspaceId: '123', orgName: 'TestOrg', workspaceName: 'TestWorkspace'],
            [workspaceId: '456', orgName: 'OtherOrg', workspaceName: 'OtherWorkspace']
        ]

        when:
        def result = cmd.getCurrentWorkspaceName(workspaces, '123')

        then:
        result == 'TestOrg / TestWorkspace'
    }

    def 'should get default workspace name when workspace does not exist'() {
        given:
        def cmd = new AuthCommandImpl()
        def workspaces = [
            [workspaceId: '123', orgName: 'TestOrg', workspaceName: 'TestWorkspace']
        ]

        when:
        def result = cmd.getCurrentWorkspaceName(workspaces, '999')

        then:
        result == 'None (Personal workspace)'
    }

    def 'should get config value from login file'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def authFile = tempDir.resolve('seqera_auth.config')
        cmd.getAuthFile() >> authFile

        def config = [:]
        def auth =['tower.accessToken': 'token-from-login']

        when:
        def result = cmd.getConfigValue(config, auth,'tower.accessToken', 'TOWER_ACCESS_TOKEN', null)

        then:
        result.value == 'token-from-login'
        result.source.endsWith('seqera_auth.config')
        result.fromConfig == true
        result.fromEnv == false
        result.isDefault == false
    }

    def 'should get config value from main config file'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def configFile = tempDir.resolve('config')
        cmd.getConfigFile() >> configFile

        def config = ['tower.endpoint': 'https://example.com']
        def auth =[:]

        when:
        def result = cmd.getConfigValue(config, auth,'tower.endpoint', 'TOWER_API_ENDPOINT', null)

        then:
        result.value == 'https://example.com'
        result.source.endsWith('config')
        result.fromConfig == true
        result.fromEnv == false
        result.isDefault == false
    }

    def 'should get config value from environment variable'() {
        given:
        def cmd = Spy(AuthCommandImpl)

        def config = [:]
        def auth =[:]

        and:
        SysEnv.push( ['TOWER_ACCESS_TOKEN': 'token-from-env'])

        when:
        def result = cmd.getConfigValue(config, auth,'tower.accessToken', 'TOWER_ACCESS_TOKEN', null)


        then:
        result.value == 'token-from-env'
        result.source == 'env var $TOWER_ACCESS_TOKEN'
        result.fromConfig == false
        result.fromEnv == true
        result.isDefault == false

        cleanup:
        SysEnv.pop()
    }

    def 'should get config value from default'() {
        given:
        def cmd = new AuthCommandImpl()

        def config = [:]
        def auth =[:]

        when:
        def result = cmd.getConfigValue(config, auth,'tower.endpoint', null, 'https://default.example.com')

        then:
        result.value == 'https://default.example.com'
        result.source == 'default'
        result.fromConfig == false
        result.fromEnv == false
        result.isDefault == true
    }

    def 'should prioritize login over config and env'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def authFile = tempDir.resolve('seqera_auth.config')
        cmd.getAuthFile() >> authFile

        def config = ['tower.accessToken': 'token-from-config']
        def auth =['tower.accessToken': 'token-from-login']

        SysEnv.push(['TOWER_ACCESS_TOKEN': 'token-from-env'])

        when:
        def result = cmd.getConfigValue(config, auth,'tower.accessToken', 'TOWER_ACCESS_TOKEN', 'default-token')

        then:
        result.value == 'token-from-login'
        result.source.endsWith('seqera_auth.config')

        cleanup:
        SysEnv.pop()
    }

    def 'should prioritize config over env when login is empty'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def configFile = tempDir.resolve('config')
        cmd.getConfigFile() >> configFile

        def config = ['tower.endpoint': 'https://config.example.com']
        def auth =[:]

        SysEnv.push(['TOWER_API_ENDPOINT': 'https://env.example.com'])

        when:
        def result = cmd.getConfigValue(config, auth,'tower.endpoint', 'TOWER_API_ENDPOINT', 'https://default.example.com')


        then:
        result.value == 'https://config.example.com'
        result.source.endsWith('config')

        cleanup:
        SysEnv.pop()
    }

    def 'should handle null environment variable name'() {
        given:
        def cmd = new AuthCommandImpl()

        def config = ['tower.enabled': true]
        def auth =[:]

        when:
        def result = cmd.getConfigValue(config, auth,'tower.enabled', null, 'false')

        then:
        result.value == true
        result.fromConfig == true
        result.fromEnv == false
    }

    def 'should print status table with simple values'() {
        given:
        def cmd = new AuthCommandImpl()
        def rows = [
            ['Setting1', 'Value1', 'Source1'],
            ['Setting2', 'Value2', 'Source2'],
            ['Setting3', 'Value3', 'Source3']
        ]

        when:
        cmd.printStatusTable(rows)
        def output = capture.toString()

        then:
        output.contains('Setting')
        output.contains('Value')
        output.contains('Source')
        output.contains('Setting1')
        output.contains('Value1')
        output.contains('Source1')
        output.contains('Setting2')
        output.contains('Value2')
        output.contains('Source2')
        output.contains('Setting3')
        output.contains('Value3')
        output.contains('Source3')
        output.contains('---') // Table separator
    }

    def 'should print status table with colored values'() {
        given:
        def cmd = new AuthCommandImpl()
        def rows = [
            ['API endpoint', nextflow.cli.ColorUtil.colorize('https://api.cloud.seqera.io', 'magenta'), 'config'],
            ['Authentication', nextflow.cli.ColorUtil.colorize('OK', 'green'), 'env var']
        ]

        when:
        cmd.printStatusTable(rows)
        def output = capture.toString()

        then:
        output.contains('API endpoint')
        output.contains('https://api.cloud.seqera.io')
        output.contains('Authentication')
        output.contains('OK')
        output.contains('config')
        output.contains('env var')
    }

    def 'should handle empty status table'() {
        given:
        def cmd = new AuthCommandImpl()
        def rows = []

        when:
        cmd.printStatusTable(rows)
        def output = capture.toString()

        then:
        output.isEmpty()
    }

    def 'should handle null status table'() {
        given:
        def cmd = new AuthCommandImpl()

        when:
        cmd.printStatusTable(null)
        def output = capture.toString()

        then:
        output.isEmpty()
    }

    def 'should strip ANSI codes correctly'() {
        given:
        def cmd = new AuthCommandImpl()
        def textWithAnsi = '\u001B[31mRed Text\u001B[0m'
        def textWithMultipleAnsi = '\u001B[1;32mBold Green\u001B[0m \u001B[34mBlue\u001B[0m'

        when:
        def stripped1 = cmd.stripAnsiCodes(textWithAnsi)
        def stripped2 = cmd.stripAnsiCodes(textWithMultipleAnsi)

        then:
        stripped1 == 'Red Text'
        stripped2 == 'Bold Green Blue'
    }

    def 'should handle null text when stripping ANSI codes'() {
        given:
        def cmd = new AuthCommandImpl()

        when:
        def result = cmd.stripAnsiCodes(null)

        then:
        result == ''
    }

    def 'should handle empty text when stripping ANSI codes'() {
        given:
        def cmd = new AuthCommandImpl()

        when:
        def result = cmd.stripAnsiCodes('')

        then:
        result == ''
    }

    def 'should pad string with ANSI codes correctly'() {
        given:
        def cmd = new AuthCommandImpl()
        def plainText = 'Hello'
        def coloredText = '\u001B[32mHello\u001B[0m'

        when:
        def paddedPlain = cmd.padStringWithAnsi(plainText, 10)
        def paddedColored = cmd.padStringWithAnsi(coloredText, 10)

        then:
        // Plain text should be padded to 10 characters
        cmd.stripAnsiCodes(paddedPlain).length() == 10
        // Colored text should preserve ANSI codes but pad visible text to 10 characters
        cmd.stripAnsiCodes(paddedColored).length() == 10
        paddedColored.contains('\u001B[32m') // ANSI codes preserved
    }

    def 'should not pad if text is already at target width'() {
        given:
        def cmd = new AuthCommandImpl()
        def text = 'Hello'

        when:
        def result = cmd.padStringWithAnsi(text, 5)

        then:
        result == 'Hello'
    }

    def 'should not pad if text exceeds target width'() {
        given:
        def cmd = new AuthCommandImpl()
        def text = 'Hello World'

        when:
        def result = cmd.padStringWithAnsi(text, 5)

        then:
        result == 'Hello World'
    }

    def 'should shorten path with user home'() {
        given:
        def cmd = new AuthCommandImpl()
        def userHome = System.getProperty('user.home')
        def path = "${userHome}/some/path/config"

        when:
        def result = cmd.shortenPath(path)

        then:
        result == '~/some/path/config'
    }

    def 'should not shorten path without user home'() {
        given:
        def cmd = new AuthCommandImpl()
        def path = '/etc/config'

        when:
        def result = cmd.shortenPath(path)

        then:
        result == '/etc/config'
    }

    def 'should print status table with varying column widths'() {
        given:
        def cmd = new AuthCommandImpl()
        def rows = [
            ['Short', 'V', 'S'],
            ['Very Long Setting Name', 'Very Long Value Here', 'Very Long Source']
        ]

        when:
        cmd.printStatusTable(rows)
        def output = capture.toString()

        then:
        // Should handle different column widths properly
        output.contains('Short')
        output.contains('Very Long Setting Name')
        output.contains('Very Long Value Here')
        output.contains('Very Long Source')
        // All rows should be aligned
        def lines = output.split('\n')
        lines.size() >= 4 // Header + separator + 2 data rows
    }

    def 'should apply minimum column widths'() {
        given:
        def cmd = new AuthCommandImpl()
        def rows = [
            ['A', 'B', 'C']
        ]

        when:
        cmd.printStatusTable(rows)
        def output = capture.toString()

        then:
        // Even with short values, should apply minimum widths
        output.contains('Setting') // Header
        output.contains('Value')
        output.contains('Source')
        // Columns should be padded to at least minimum width
        def lines = output.split('\n')
        lines.size() >= 3 // Header + separator + data row
    }

    def 'should collect status with valid authentication'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def config = [:]
        def authConfig = [
            'tower.accessToken': 'test-token',
            'tower.endpoint': 'https://api.cloud.seqera.io',
            'tower.enabled': true
        ]

        // Mock API calls
        cmd.checkApiConnection(_) >> true
        cmd.callUserInfoApi(_, _) >> [userName: 'testuser', id: '123']
        cmd.getWorkspaceDetailsFromApi(_, _, _) >> null

        when:
        def status = cmd.collectStatus(config, authConfig)

        then:
        status != null
        status.table.size() == 5 // endpoint, connection, auth, monitoring, workspace
        // Check API endpoint row
        status.table[0][0] == 'API endpoint'
        status.table[0][1].contains('https://api.cloud.seqera.io')
        // Check API connection row
        status.table[1][0] == 'API connection'
        status.table[1][1].contains('OK')
        // Check authentication row
        status.table[2][0] == 'Authentication'
        status.table[2][1].contains('OK')
        status.table[2][1].contains('testuser')
        // Check monitoring row
        status.table[3][0] == 'Workflow monitoring'
        status.table[3][1].contains('Yes')
        // Check workspace row
        status.table[4][0] == 'Default workspace'
        status.table[4][1].contains('Personal workspace')
    }

    def 'should collect status without authentication'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def config = [:]
        def authConfig = [:]

        cmd.checkApiConnection(_) >> true

        when:
        def status = cmd.collectStatus(config, authConfig)

        then:
        status != null
        status.table.size() == 5
        // Authentication should show error
        status.table[2][0] == 'Authentication'
        status.table[2][1].contains('ERROR')
        status.table[2][1].contains('no token')
        status.table[2][2] == 'not set'
    }

    def 'should collect status with failed API connection'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def config = [:]
        def authConfig = ['tower.endpoint': 'https://unreachable.example.com']

        cmd.checkApiConnection(_) >> false

        when:
        def status = cmd.collectStatus(config, authConfig)

        then:
        status != null
        status.table[1][0] == 'API connection'
        status.table[1][1].contains('ERROR')
    }

    def 'should collect status with failed authentication'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def config = [:]
        def authConfig = ['tower.accessToken': 'invalid-token']

        cmd.checkApiConnection(_) >> true
        cmd.callUserInfoApi(_, _) >> { throw new RuntimeException('Invalid token') }

        when:
        def status = cmd.collectStatus(config, authConfig)

        then:
        status != null
        status.table[2][0] == 'Authentication'
        status.table[2][1].contains('ERROR')
        status.table[2][2] == 'failed'
    }

    def 'should collect status with monitoring enabled'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def config = [:]
        def authConfig = ['tower.enabled': true]

        cmd.checkApiConnection(_) >> true

        when:
        def status = cmd.collectStatus(config, authConfig)

        then:
        status != null
        status.table[3][0] == 'Workflow monitoring'
        status.table[3][1].contains('Yes')
    }

    def 'should collect status with monitoring disabled'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def config = [:]
        def authConfig = ['tower.enabled': false]

        cmd.checkApiConnection(_) >> true

        when:
        def status = cmd.collectStatus(config, authConfig)

        then:
        status != null
        status.table[3][0] == 'Workflow monitoring'
        status.table[3][1].contains('No')
    }

    def 'should collect status with workspace details'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def config = [:]
        def authConfig = [
            'tower.accessToken': 'test-token',
            'tower.workspaceId': '12345'
        ]

        cmd.checkApiConnection(_) >> true
        cmd.callUserInfoApi(_, _) >> [userName: 'testuser', id: '123']
        cmd.getWorkspaceDetailsFromApi(_, _, _) >> [
            orgName: 'TestOrg',
            workspaceName: 'TestWorkspace',
            workspaceFullName: 'test-org/test-workspace'
        ]

        when:
        def status = cmd.collectStatus(config, authConfig)

        then:
        status != null
        status.table[4][0] == 'Default workspace'
        status.table[4][1].contains('12345')
        status.workspaceInfo != null
        status.workspaceInfo.orgName == 'TestOrg'
        status.workspaceInfo.workspaceName == 'TestWorkspace'
    }

    def 'should collect status with workspace ID but no details'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def config = [:]
        def authConfig = [
            'tower.accessToken': 'test-token',
            'tower.workspaceId': '12345'
        ]

        cmd.checkApiConnection(_) >> true
        cmd.callUserInfoApi(_, _) >> [userName: 'testuser', id: '123']
        cmd.getWorkspaceDetailsFromApi(_, _, _) >> null

        when:
        def status = cmd.collectStatus(config, authConfig)

        then:
        status != null
        status.table[4][0] == 'Default workspace'
        status.table[4][1].contains('12345')
        status.workspaceInfo == null
    }

    def 'should collect status from environment variables'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def config = [:]
        def authConfig = [:]

        cmd.checkApiConnection(_) >> true
        cmd.callUserInfoApi(_, _) >> [userName: 'envuser', id: '456']
        cmd.getWorkspaceDetailsFromApi(_,_,_) >> [:]

        SysEnv.push(['TOWER_ACCESS_TOKEN': 'env-token',
                     'TOWER_API_ENDPOINT': 'https://env.example.com',
                     'TOWER_WORKFLOW_ID': 'ws-123'] )

        when:
        def status = cmd.collectStatus(config, authConfig)


        then:
        status != null
        status.table[0][1].contains('https://env.example.com')
        status.table[0][2] == 'env var $TOWER_API_ENDPOINT'
        status.table[2][1].contains('envuser')
        status.table[2][2].contains('env var $TOWER_ACCESS_TOKEN')
        status.table[4][2].contains('env var $TOWER_WORKFLOW_ID')

        cleanup:
        SysEnv.pop()
    }

    def 'should collect status with default values'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def config = [:]
        def authConfig = [:]

        cmd.checkApiConnection(_) >> true

        when:
        def status = cmd.collectStatus(config, authConfig)

        then:
        status != null
        // Should use default endpoint
        status.table[0][1].contains('https://api.cloud.seqera.io')
        status.table[0][2] == 'default'
        // Should show monitoring disabled by default
        status.table[3][1].contains('No')
        status.table[3][2] == 'default'
    }

    def 'should collect status with mixed sources'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def authFile = tempDir.resolve('seqera_auth.config')
        def configFile = tempDir.resolve('config')

        cmd.getAuthFile() >> authFile
        cmd.getConfigFile() >> configFile

        def config = ['tower.enabled': true]
        def authConfig = ['tower.accessToken': 'login-token']

        cmd.checkApiConnection(_) >> true
        cmd.callUserInfoApi(_, _) >> [userName: 'mixeduser', id: '789']
        SysEnv.push(['TOWER_WORKFLOW_ID': 'ws-env'])

        when:
        def status = cmd.collectStatus(config, authConfig)


        then:
        status != null
        // Token from login file
        status.table[2][2].endsWith('.login')
        // Enabled from config file
        status.table[3][2].endsWith('config')
        // Workspace from env var
        status.table[4][2].contains('env var $TOWER_WORKFLOW_ID')

        cleanup:
        SysEnv.pop()
    }

    def 'should print status correctly'() {
        given:
        def cmd = new AuthCommandImpl()
        def status = new AuthCommandImpl.ConfigStatus(
            [
                ['API endpoint', 'https://api.cloud.seqera.io', 'default'],
                ['Authentication', 'OK', 'config']
            ],
            [
                orgName: 'TestOrg',
                workspaceName: 'TestWorkspace',
                workspaceFullName: 'test-org/test-workspace'
            ]
        )

        when:
        cmd.printStatus(status)
        def output = capture.toString()

        then:
        output.contains('API endpoint')
        output.contains('https://api.cloud.seqera.io')
        output.contains('Authentication')
        output.contains('OK')
        output.contains('TestOrg')
        output.contains('TestWorkspace')
    }

    def 'should print status without workspace info'() {
        given:
        def cmd = new AuthCommandImpl()
        def status = new AuthCommandImpl.ConfigStatus(
            [
                ['API endpoint', 'https://api.cloud.seqera.io', 'default']
            ],
            null
        )

        when:
        cmd.printStatus(status)
        def output = capture.toString()

        then:
        output.contains('API endpoint')
        !output.contains('TestOrg')
    }

    def 'should detect existing auth file and prevent duplicate login'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def authFile = tempDir.resolve('seqera_auth.config')
        Files.createFile(authFile)

        cmd.getAuthFile() >> authFile

        when:
        cmd.login('https://api.cloud.seqera.io')
        def output = capture.toString()

        then:
        output.contains('Error: Authentication token is already configured')
        output.contains('nextflow auth logout')
    }

    def 'should warn when TOWER_ACCESS_TOKEN env var is set during login'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def authFile = tempDir.resolve('seqera_auth.config-not-exists')

        cmd.getAuthFile() >> authFile
        cmd.performAuth0Login(_, _) >> { /* mock to prevent actual login */ }
        SysEnv.push(['TOWER_ACCESS_TOKEN': 'env-token'])

        when:
        cmd.login('https://api.cloud.seqera.io')
        def output = capture.toString()

        then:
        output.contains('WARNING: Authentication token is already configured via TOWER_ACCESS_TOKEN environment variable')

        cleanup:
        SysEnv.pop()
    }

    def 'should normalize API URL during login'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def authFile = tempDir.resolve('seqera_auth.config-not-exists')

        cmd.getAuthFile() >> authFile

        when:
        cmd.login('api.cloud.seqera.io')

        then:
        1 * cmd.performAuth0Login('https://api.cloud.seqera.io', _) >> { /* mock to prevent actual login */ }
    }

    def 'should route to enterprise auth for non-cloud endpoints'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def authFile = tempDir.resolve('seqera_auth.config-not-exists')

        cmd.getAuthFile() >> authFile

        when:
        cmd.login('https://enterprise.example.com')

        then:
        1 * cmd.handleEnterpriseAuth('https://enterprise.example.com') >> {/* mock to prevent actual login */ }
        0 * cmd.performAuth0Login(_, _)
    }

    def 'should route to Auth0 login for cloud endpoints'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def authFile = tempDir.resolve('seqera_auth.config-not-exists')

        cmd.getAuthFile() >> authFile

        when:
        cmd.login('https://api.cloud.seqera.io')

        then:
        1 * cmd.performAuth0Login('https://api.cloud.seqera.io', _) >>  { /* mock to prevent actual login */ }
        0 * cmd.handleEnterpriseAuth(_)
    }

    def 'should save auth to config after successful PAT generation'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def authFile = tempDir.resolve('seqera_auth.config')
        def configFile = tempDir.resolve('config')

        cmd.getAuthFile() >> authFile
        cmd.getConfigFile() >> configFile

        def config = [
            'tower.accessToken': 'generated-pat-token',
            'tower.endpoint': 'https://api.cloud.seqera.io',
            'tower.enabled': true
        ]

        when:
        cmd.saveAuthToConfig('generated-pat-token', 'https://api.cloud.seqera.io')

        then:
        Files.exists(authFile)
        def content = Files.readString(authFile)
        content.contains('accessToken = \'generated-pat-token\'')
        content.contains('endpoint = \'https://api.cloud.seqera.io\'')
        content.contains('enabled = true')
    }

    def 'should perform Auth0 request correctly'() {
        given:
        def cmd = Spy(AuthCommandImpl)

        cmd.createHttpClient(_) >> {
            // Mock HxClient that returns successful response
            def mockClient = Mock(HxClient)
            def mockResponse = Mock(HttpResponse)
            mockResponse.statusCode() >> 200
            mockResponse.body() >> '{"access_token":"test-token","expires_in":3600}'
            mockClient.send(_, _) >> mockResponse
            return mockClient
        }

        def params = [
            'client_id': 'test-client-id',
            'scope': 'openid profile'
        ]

        when:
        def result = cmd.performAuth0Request('https://auth.example.com/oauth/token', params)

        then:
        result != null
        result['access_token'] == 'test-token'
        result['expires_in'] == 3600
    }

    def 'should handle Auth0 request errors'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        cmd.createHttpClient(_) >> {
            def mockClient = Mock(HxClient)
            def mockResponse = Mock(HttpResponse)
            mockResponse.statusCode() >> 400
            mockResponse.body() >> '{"error":"invalid_grant","error_description":"Invalid credentials"}'
            mockClient.send(_, _) >> mockResponse
            return mockClient
        }

        def params = ['client_id': 'test-client-id']

        when:
        cmd.performAuth0Request('https://auth.example.com/oauth/token', params)

        then:
        def ex = thrown(RuntimeException)
        ex.message.contains('invalid_grant')
        ex.message.contains('Invalid credentials')
    }

    def 'should generate PAT with correct token name'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def username = System.getProperty('user.name')

        cmd.createHttpClient(_) >> {
            def mockClient = Mock(HxClient)
            def mockResponse = Mock(HttpResponse)
            mockResponse.statusCode() >> 200
            mockResponse.body() >> '{"accessKey":"generated-pat-123","id":"token-id-456"}'
            mockClient.send(_, _) >> mockResponse
            return mockClient
        }

        when:
        def pat = cmd.generatePAT('auth-token', 'https://api.cloud.seqera.io')

        then:
        pat == 'generated-pat-123'
    }

    def 'should fail PAT generation on error response'() {
        given:
        def cmd = Spy(AuthCommandImpl)

        cmd.createHttpClient(_) >> {
            def mockClient = Mock(HxClient)
            def mockResponse = Mock(HttpResponse)
            mockResponse.statusCode() >> 401
            mockResponse.body() >> 'Unauthorized'
            mockClient.send(_, _) >> mockResponse
            return mockClient
        }

        when:
        cmd.generatePAT('invalid-token', 'https://api.cloud.seqera.io')

        then:
        def ex = thrown(RuntimeException)
        ex.message.contains('Failed to generate PAT')
    }

    def 'should request device authorization with correct parameters'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def auth0Config = [
            domain: 'seqera.eu.auth0.com',
            clientId: 'test-client-id'
        ]

        cmd.performAuth0Request(_, _) >> [
            device_code: 'device-code-123',
            user_code: 'ABCD-1234',
            verification_uri: 'https://seqera.eu.auth0.com/activate',
            interval: 5
        ]

        when:
        def result = cmd.requestDeviceAuthorization(auth0Config)

        then:
        result['device_code'] == 'device-code-123'
        result['user_code'] == 'ABCD-1234'
        result['verification_uri'] == 'https://seqera.eu.auth0.com/activate'
    }

    def 'should poll for device token and return on success'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def auth0Config = [
            domain: 'seqera.eu.auth0.com',
            clientId: 'test-client-id'
        ]

        // First call returns pending, second call returns token
        int callCount = 0
        cmd.performAuth0Request(_, _) >> {
            if (callCount++ == 0) {
                throw new RuntimeException('authorization_pending')
            } else {
                return [access_token: 'final-token', token_type: 'Bearer']
            }
        }

        when:
        def result = cmd.pollForDeviceToken('device-code-123', 1, auth0Config)

        then:
        result['access_token'] == 'final-token'
    }

    def 'should handle expired token during polling'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def auth0Config = [
            domain: 'seqera.eu.auth0.com',
            clientId: 'test-client-id'
        ]

        cmd.performAuth0Request(_, _) >> {
            throw new RuntimeException('expired_token')
        }

        when:
        cmd.pollForDeviceToken('device-code-123', 1, auth0Config)

        then:
        def ex = thrown(RuntimeException)
        ex.message.contains('expired')
    }

    def 'should handle access denied during polling'() {
        given:
        def cmd = Spy(AuthCommandImpl)
        def auth0Config = [
            domain: 'seqera.eu.auth0.com',
            clientId: 'test-client-id'
        ]

        cmd.performAuth0Request(_, _) >> {
            throw new RuntimeException('access_denied')
        }

        when:
        cmd.pollForDeviceToken('device-code-123', 1, auth0Config)

        then:
        def ex = thrown(RuntimeException)
        ex.message.contains('denied')
    }
}
