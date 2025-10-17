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

package io.seqera.tower.plugin.launch

import nextflow.cli.CmdLaunch
import nextflow.exception.AbortOperationException
import org.junit.Rule
import spock.lang.Specification
import spock.lang.TempDir
import test.OutputCapture

import java.nio.file.Files
import java.nio.file.Path

/**
 * Test LaunchCommandImpl functionality
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class LaunchCommandImplTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @TempDir
    Path tempDir

    // ===== Pipeline Validation Tests =====

    def 'should reject local file path'() {
        given:
        def cmd = new LaunchCommandImpl()

        when:
        cmd.validateAndResolvePipeline('/path/to/local/workflow')

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('Local file paths are not supported')
    }

    def 'should reject relative file path with ./'() {
        given:
        def cmd = new LaunchCommandImpl()

        when:
        cmd.validateAndResolvePipeline('./local/workflow')

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('Local file paths are not supported')
    }

    def 'should reject relative file path with ../'() {
        given:
        def cmd = new LaunchCommandImpl()

        when:
        cmd.validateAndResolvePipeline('../local/workflow')

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('Local file paths are not supported')
    }

    def 'should accept remote repository URL'() {
        given:
        def cmd = Spy(LaunchCommandImpl)

        and:
        cmd.resolvePipelineUrl(_) >> 'https://github.com/org/repo'

        when:
        def result = cmd.validateAndResolvePipeline('https://github.com/org/repo')

        then:
        result == 'https://github.com/org/repo'
    }

    def 'should accept github short name'() {
        given:
        def cmd = Spy(LaunchCommandImpl)

        and:
        cmd.resolvePipelineUrl(_) >> 'https://github.com/nf-core/rnaseq'

        when:
        def result = cmd.validateAndResolvePipeline('nf-core/rnaseq')

        then:
        result == 'https://github.com/nf-core/rnaseq'
    }

    def 'should identify local path correctly'() {
        given:
        def cmd = new LaunchCommandImpl()

        expect:
        cmd.isLocalPath('/absolute/path') == true
        cmd.isLocalPath('./relative/path') == true
        cmd.isLocalPath('../parent/path') == true
        cmd.isLocalPath('C:\\windows\\path') == true
        cmd.isLocalPath('remote/repo') == false
        cmd.isLocalPath('https://github.com/org/repo') == false
    }

    // ===== Parameter Parsing Tests =====

    def 'should parse simple parameters'() {
        given:
        def cmd = new LaunchCommandImpl()
        def params = ['input': 'data.csv', 'output': 'results/']

        when:
        def paramsText = cmd.buildParamsText(params, null)

        then:
        paramsText != null
        paramsText.contains('"input"')
        paramsText.contains('data.csv')
        paramsText.contains('"output"')
        paramsText.contains('results/')
    }

    def 'should parse nested parameters with dot notation'() {
        given:
        def cmd = new LaunchCommandImpl()
        def params = ['genome.fasta': 'hg38.fa', 'genome.index': 'hg38.idx']

        when:
        def paramsText = cmd.buildParamsText(params, null)

        then:
        paramsText != null
        paramsText.contains('"genome"')
        paramsText.contains('"fasta"')
        paramsText.contains('hg38.fa')
        paramsText.contains('"index"')
        paramsText.contains('hg38.idx')
    }

    def 'should parse boolean parameters'() {
        given:
        def cmd = new LaunchCommandImpl()
        def params = ['verbose': 'true', 'skip': 'false']

        when:
        def paramsText = cmd.buildParamsText(params, null)

        then:
        paramsText != null
        paramsText.contains('"verbose":true')
        paramsText.contains('"skip":false')
    }

    def 'should parse numeric parameters'() {
        given:
        def cmd = new LaunchCommandImpl()
        def params = ['threads': '8', 'memory': '16.5', 'size': '1000000']

        when:
        def paramsText = cmd.buildParamsText(params, null)

        then:
        paramsText != null
        paramsText.contains('"threads":8')
        paramsText.contains('"memory":16.5')
        paramsText.contains('"size":1000000')
    }

    def 'should parse parameters from JSON file'() {
        given:
        def cmd = new LaunchCommandImpl()
        def paramsFile = tempDir.resolve('params.json')
        Files.writeString(paramsFile, '{"input": "data.csv", "output": "results/"}')

        when:
        def paramsText = cmd.buildParamsText([:], paramsFile.toString())

        then:
        paramsText != null
        paramsText.contains('"input"')
        paramsText.contains('data.csv')
    }

    def 'should parse parameters from YAML file'() {
        given:
        def cmd = new LaunchCommandImpl()
        def paramsFile = tempDir.resolve('params.yml')
        Files.writeString(paramsFile, 'input: data.csv\noutput: results/')

        when:
        def paramsText = cmd.buildParamsText([:], paramsFile.toString())

        then:
        paramsText != null
        paramsText.contains('"input"')
        paramsText.contains('data.csv')
    }

    def 'should merge CLI params with params file, CLI taking precedence'() {
        given:
        def cmd = new LaunchCommandImpl()
        def paramsFile = tempDir.resolve('params.json')
        Files.writeString(paramsFile, '{"input": "file.csv", "output": "old/"}')
        def params = ['output': 'new/']

        when:
        def paramsText = cmd.buildParamsText(params, paramsFile.toString())

        then:
        paramsText != null
        paramsText.contains('"input"')
        paramsText.contains('file.csv')
        paramsText.contains('"output"')
        paramsText.contains('new/')
    }

    def 'should reject invalid params file extension'() {
        given:
        def cmd = new LaunchCommandImpl()
        def paramsFile = tempDir.resolve('params.txt')
        Files.writeString(paramsFile, 'input: data.csv')

        when:
        cmd.buildParamsText([:], paramsFile.toString())

        then:
        thrown(AbortOperationException)
    }

    def 'should handle missing params file'() {
        given:
        def cmd = new LaunchCommandImpl()

        when:
        cmd.buildParamsText([:], '/nonexistent/params.json')

        then:
        thrown(AbortOperationException)
    }

    def 'should return null when no parameters provided'() {
        given:
        def cmd = new LaunchCommandImpl()

        when:
        def paramsText = cmd.buildParamsText([:], null)

        then:
        paramsText == null
    }

    def 'should convert kebab-case to camelCase in parameter names'() {
        given:
        def cmd = new LaunchCommandImpl()
        def params = ['max-memory': '16GB', 'output-dir': 'results/']

        when:
        def paramsText = cmd.buildParamsText(params, null)

        then:
        paramsText != null
        paramsText.contains('"maxMemory"')
        paramsText.contains('"outputDir"')
    }

    def 'should handle escaped dots in parameter names'() {
        given:
        def cmd = new LaunchCommandImpl()
        def params = ['file\\.name': 'test.csv']

        when:
        def paramsText = cmd.buildParamsText(params, null)

        then:
        paramsText != null
        paramsText.contains('"file.name"')
    }

    // ===== Config File Tests =====

    def 'should read config file content'() {
        given:
        def cmd = new LaunchCommandImpl()
        def configFile = tempDir.resolve('test.config')
        Files.writeString(configFile, 'process.cpus = 8\nprocess.memory = "16 GB"')

        when:
        def configText = cmd.buildConfigText([configFile.toString()])

        then:
        configText != null
        configText.contains('process.cpus = 8')
        configText.contains('process.memory = "16 GB"')
    }

    def 'should concatenate multiple config files'() {
        given:
        def cmd = new LaunchCommandImpl()
        def config1 = tempDir.resolve('config1.config')
        def config2 = tempDir.resolve('config2.config')
        Files.writeString(config1, 'process.cpus = 8')
        Files.writeString(config2, 'process.memory = "16 GB"')

        when:
        def configText = cmd.buildConfigText([config1.toString(), config2.toString()])

        then:
        configText != null
        configText.contains('process.cpus = 8')
        configText.contains('process.memory = "16 GB"')
    }

    def 'should return null when no config files provided'() {
        given:
        def cmd = new LaunchCommandImpl()

        when:
        def configText = cmd.buildConfigText(null)

        then:
        configText == null

        when:
        configText = cmd.buildConfigText([])

        then:
        configText == null
    }

    def 'should throw exception for missing config file'() {
        given:
        def cmd = new LaunchCommandImpl()

        when:
        cmd.buildConfigText(['/nonexistent/config.config'])

        then:
        thrown(AbortOperationException)
    }

    // ===== Launch Context Tests =====

    def 'should throw error when no access token configured'() {
        given:
        def cmd = Spy(LaunchCommandImpl)
        cmd.readConfig() >> [:]
        def options = new CmdLaunch.LaunchOptions(pipeline: 'nf-core/rnaseq')

        when:
        cmd.initializeLaunchContext(options)

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('No authentication found')
        ex.message.contains('nextflow auth login')
    }

    def 'should initialize context with valid config'() {
        given:
        def cmd = Spy(LaunchCommandImpl)
        def config = [
            'tower.accessToken': 'test-token',
            'tower.endpoint': 'https://api.cloud.seqera.io'
        ]
        cmd.readConfig() >> config
        cmd.getUserInfo(_, _) >> [name: 'testuser', id: '123']
        cmd.resolveWorkspaceId(_, _, _, _) >> null
        cmd.getWorkspaceDetails(_, _, _) >> null
        cmd.resolveComputeEnvironment(_, _, _, _) >> [id: 'ce-123', name: 'test-ce', workDir: 's3://bucket/work']

        def options = new CmdLaunch.LaunchOptions(pipeline: 'nf-core/rnaseq')

        when:
        def context = cmd.initializeLaunchContext(options)

        then:
        context.accessToken == 'test-token'
        context.apiEndpoint == 'https://api.cloud.seqera.io'
        context.userName == 'testuser'
        context.computeEnvId == 'ce-123'
        context.computeEnvName == 'test-ce'
        context.workDir == 's3://bucket/work'
    }

    def 'should use default endpoint when not configured'() {
        given:
        def cmd = Spy(LaunchCommandImpl)
        def config = ['tower.accessToken': 'test-token']
        cmd.readConfig() >> config
        cmd.getUserInfo(_, _) >> [name: 'testuser', id: '123']
        cmd.resolveWorkspaceId(_, _, _, _) >> null
        cmd.getWorkspaceDetails(_, _, _) >> null
        cmd.resolveComputeEnvironment(_, _, _, _) >> [id: 'ce-123', name: 'test-ce', workDir: 's3://bucket/work']

        def options = new CmdLaunch.LaunchOptions(pipeline: 'nf-core/rnaseq')

        when:
        def context = cmd.initializeLaunchContext(options)

        then:
        context.apiEndpoint == 'https://api.cloud.seqera.io'
    }

    def 'should resolve workspace details when workspace ID provided'() {
        given:
        def cmd = Spy(LaunchCommandImpl)
        def config = ['tower.accessToken': 'test-token', 'tower.workspaceId': 12345]
        cmd.readConfig() >> config
        cmd.getUserInfo(_, _) >> [name: 'testuser', id: '123']
        cmd.resolveWorkspaceId(_, _, _, _) >> 12345L
        cmd.getWorkspaceDetails(_, _, _) >> [orgName: 'TestOrg', workspaceName: 'TestWS']
        cmd.resolveComputeEnvironment(_, _, _, _) >> [id: 'ce-123', name: 'test-ce', workDir: 's3://bucket/work']

        def options = new CmdLaunch.LaunchOptions(pipeline: 'nf-core/rnaseq')

        when:
        def context = cmd.initializeLaunchContext(options)

        then:
        context.workspaceId == 12345L
        context.orgName == 'TestOrg'
        context.workspaceName == 'TestWS'
    }

    // ===== Compute Environment Tests =====

    def 'should find compute environment by name'() {
        given:
        def cmd = Spy(LaunchCommandImpl)
        def computeEnvs = [
            [id: 'ce-1', name: 'primary-ce', primary: true],
            [id: 'ce-2', name: 'secondary-ce', primary: false]
        ]
        cmd.getComputeEnvironments(_, _, _) >> computeEnvs

        when:
        def result = cmd.findComputeEnv('secondary-ce', null, 'token', 'endpoint')

        then:
        result.id == 'ce-2'
        result.name == 'secondary-ce'
    }

    def 'should find primary compute environment when name not provided'() {
        given:
        def cmd = Spy(LaunchCommandImpl)
        def computeEnvs = [
            [id: 'ce-1', name: 'primary-ce', primary: true],
            [id: 'ce-2', name: 'secondary-ce', primary: false]
        ]
        cmd.getComputeEnvironments(_, _, _) >> computeEnvs

        when:
        def result = cmd.findComputeEnv(null, null, 'token', 'endpoint')

        then:
        result.id == 'ce-1'
        result.name == 'primary-ce'
        result.primary == true
    }

    def 'should return null when compute environment not found'() {
        given:
        def cmd = Spy(LaunchCommandImpl)
        cmd.getComputeEnvironments(_, _, _) >> []

        when:
        def result = cmd.findComputeEnv('nonexistent', null, 'token', 'endpoint')

        then:
        result == null
    }

    def 'should throw error when compute environment not found'() {
        given:
        def cmd = Spy(LaunchCommandImpl)
        // Mock findComputeEnv to return null (not found)
        cmd.findComputeEnv('nonexistent', null, 'token', 'https://api.cloud.seqera.io') >> null

        when:
        cmd.resolveComputeEnvironment('nonexistent', null, 'token', 'https://api.cloud.seqera.io')

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('Compute environment \'nonexistent\' not found')
    }

    def 'should throw error when no primary compute environment'() {
        given:
        def cmd = Spy(LaunchCommandImpl)
        // Mock findComputeEnv to return null (no primary found)
        cmd.findComputeEnv(null, null, 'token', 'https://api.cloud.seqera.io') >> null

        when:
        cmd.resolveComputeEnvironment(null, null, 'token', 'https://api.cloud.seqera.io')

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('No primary compute environment found')
    }

    // ===== Work Directory Tests =====

    def 'should use CLI work dir when provided'() {
        given:
        def cmd = new LaunchCommandImpl()
        def computeEnvInfo = [workDir: 's3://default/work']

        when:
        def workDir = cmd.resolveWorkDirectory('s3://custom/work', computeEnvInfo)

        then:
        workDir == 's3://custom/work'
    }

    def 'should use compute environment work dir when CLI not provided'() {
        given:
        def cmd = new LaunchCommandImpl()
        def computeEnvInfo = [workDir: 's3://default/work']

        when:
        def workDir = cmd.resolveWorkDirectory(null, computeEnvInfo)

        then:
        workDir == 's3://default/work'
    }

    def 'should throw error when no work dir available'() {
        given:
        def cmd = new LaunchCommandImpl()
        def computeEnvInfo = [:]

        when:
        cmd.resolveWorkDirectory(null, computeEnvInfo)

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('Work directory is required')
    }

    // ===== Launch Request Building Tests =====

    def 'should build basic launch request'() {
        given:
        def cmd = new LaunchCommandImpl()
        def options = new CmdLaunch.LaunchOptions(pipeline: 'nf-core/rnaseq')
        def context = new LaunchCommandImpl.LaunchContext(
            computeEnvId: 'ce-123',
            workDir: 's3://bucket/work'
        )

        when:
        def request = cmd.buildLaunchRequestPayload(options, context, 'https://github.com/nf-core/rnaseq', null, null)

        then:
        request.launch.computeEnvId == 'ce-123'
        request.launch.workDir == 's3://bucket/work'
        request.launch.pipeline == 'https://github.com/nf-core/rnaseq'
        request.launch.resume == false
        request.launch.pullLatest == false
        request.launch.stubRun == false
    }

    def 'should include optional parameters in launch request'() {
        given:
        def cmd = new LaunchCommandImpl()
        def options = new CmdLaunch.LaunchOptions(
            pipeline: 'nf-core/rnaseq',
            runName: 'test-run',
            revision: 'main',
            profile: 'test',
            resume: 'session-id',
            latest: true,
            stubRun: true,
            mainScript: 'main.nf',
            entryName: 'workflow1'
        )
        def context = new LaunchCommandImpl.LaunchContext(
            computeEnvId: 'ce-123',
            workDir: 's3://bucket/work'
        )

        when:
        def request = cmd.buildLaunchRequestPayload(options, context, 'https://github.com/nf-core/rnaseq',
                                                     '{"input":"data.csv"}', 'process.cpus = 8')

        then:
        request.launch.runName == 'test-run'
        request.launch.revision == 'main'
        request.launch.configProfiles == 'test'
        request.launch.resume == true
        request.launch.pullLatest == true
        request.launch.stubRun == true
        request.launch.mainScript == 'main.nf'
        request.launch.entryName == 'workflow1'
        request.launch.paramsText == '{"input":"data.csv"}'
        request.launch.configText == 'process.cpus = 8'
    }

    // ===== Workflow Status Tests =====

    def 'should get color for workflow status'() {
        given:
        def cmd = new LaunchCommandImpl()

        expect:
        cmd.getColorForStatus('PENDING') == 'yellow'
        cmd.getColorForStatus('SUBMITTED') == 'yellow'
        cmd.getColorForStatus('RUNNING') == 'blue'
        cmd.getColorForStatus('SUCCEEDED') == 'green'
        cmd.getColorForStatus('FAILED') == 'red'
        cmd.getColorForStatus('CANCELLED') == 'red'
        cmd.getColorForStatus('ABORTED') == 'red'
        cmd.getColorForStatus(null) == 'cyan'
    }

    def 'should get spinner mode for workflow status'() {
        given:
        def cmd = new LaunchCommandImpl()

        expect:
        cmd.getSpinnerMode('PENDING', false) == 'waiting'
        cmd.getSpinnerMode('SUBMITTED', false) == 'waiting'
        cmd.getSpinnerMode('RUNNING', false) == 'running'
        cmd.getSpinnerMode('SUCCEEDED', false) == 'succeeded'
        cmd.getSpinnerMode('FAILED', false) == 'failed'
        cmd.getSpinnerMode('CANCELLED', false) == 'failed'
        cmd.getSpinnerMode('ABORTED', false) == 'failed'
        cmd.getSpinnerMode(null, false) == 'waiting'
    }

    def 'should format workflow status'() {
        given:
        def cmd = new LaunchCommandImpl()

        when:
        def formatted = cmd.formatWorkflowStatus('RUNNING')

        then:
        formatted.contains('Workflow status:')
        formatted.contains('RUNNING')
    }

    // ===== Workspace Resolution Tests =====

    def 'should use workspace ID from config'() {
        given:
        def cmd = new LaunchCommandImpl()
        def config = ['tower.workspaceId': 12345L]

        when:
        def workspaceId = cmd.resolveWorkspaceId(config, null, 'token', 'endpoint')

        then:
        workspaceId == 12345L
    }

    def 'should lookup workspace by name'() {
        given:
        def cmd = Spy(LaunchCommandImpl)
        def config = [:]
        def workspaces = [
            [workspaceId: 111, workspaceName: 'ws1'],
            [workspaceId: 222, workspaceName: 'ws2']
        ]
        cmd.getUserInfo(_, _) >> [id: 'user-123']
        cmd.getUserWorkspaces(_, _, _) >> workspaces

        when:
        def workspaceId = cmd.resolveWorkspaceId(config, 'ws2', 'token', 'endpoint')

        then:
        workspaceId == 222
    }

    def 'should throw error when workspace not found by name'() {
        given:
        def cmd = Spy(LaunchCommandImpl)
        def config = [:]
        cmd.getUserInfo(_, _) >> [id: 'user-123']
        cmd.getUserWorkspaces(_, _, _) >> []

        when:
        cmd.resolveWorkspaceId(config, 'nonexistent', 'token', 'endpoint')

        then:
        def ex = thrown(AbortOperationException)
        ex.message.contains('Workspace \'nonexistent\' not found')
    }

    def 'should return null when no workspace specified'() {
        given:
        def cmd = new LaunchCommandImpl()
        def config = [:]

        when:
        def workspaceId = cmd.resolveWorkspaceId(config, null, 'token', 'endpoint')

        then:
        workspaceId == null
    }

    // ===== URL Building Tests =====

    def 'should build URL without query params'() {
        given:
        def cmd = new LaunchCommandImpl()

        when:
        def url = cmd.buildUrl('https://api.cloud.seqera.io', '/workflow/launch', [:])

        then:
        url == 'https://api.cloud.seqera.io/workflow/launch'
    }

    def 'should build URL with query params'() {
        given:
        def cmd = new LaunchCommandImpl()

        when:
        def url = cmd.buildUrl('https://api.cloud.seqera.io', '/workflow/launch', [workspaceId: '12345'])

        then:
        url.contains('https://api.cloud.seqera.io/workflow/launch?')
        url.contains('workspaceId=12345')
    }

    def 'should URL encode query params'() {
        given:
        def cmd = new LaunchCommandImpl()

        when:
        def url = cmd.buildUrl('https://api.cloud.seqera.io', '/workflow', [name: 'test workflow'])

        then:
        url.contains('name=test+workflow')
    }

    // ===== Launch Result Tests =====

    def 'should extract launch result with workflow details'() {
        given:
        def cmd = new LaunchCommandImpl()
        def response = [workflowId: 'wf-123']
        def workflowDetails = [
            workflow: [
                runName: 'test-run',
                commitId: 'abc123',
                revision: 'main'
            ]
        ]
        def options = new CmdLaunch.LaunchOptions(pipeline: 'nf-core/rnaseq')
        def context = new LaunchCommandImpl.LaunchContext(
            apiEndpoint: 'https://api.cloud.seqera.io',
            userName: 'testuser'
        )

        when:
        def result = cmd.extractLaunchResult(response, workflowDetails, options, 'https://github.com/nf-core/rnaseq', context)

        then:
        result.workflowId == 'wf-123'
        result.runName == 'test-run'
        result.commitId == 'abc123'
        result.revision == 'main'
        result.repository == 'https://github.com/nf-core/rnaseq'
        result.trackingUrl.contains('/user/testuser/watch/wf-123/')
    }

    def 'should extract launch result without workflow details'() {
        given:
        def cmd = new LaunchCommandImpl()
        def response = [workflowId: 'wf-123']
        def options = new CmdLaunch.LaunchOptions(pipeline: 'nf-core/rnaseq', runName: 'custom-run')
        def context = new LaunchCommandImpl.LaunchContext(
            apiEndpoint: 'https://api.cloud.seqera.io',
            userName: 'testuser'
        )

        when:
        def result = cmd.extractLaunchResult(response, null, options, 'https://github.com/nf-core/rnaseq', context)

        then:
        result.workflowId == 'wf-123'
        result.runName == 'custom-run'
        result.commitId == 'unknown'
    }

    def 'should build tracking URL for organization workspace'() {
        given:
        def cmd = new LaunchCommandImpl()
        def response = [workflowId: 'wf-123']
        def options = new CmdLaunch.LaunchOptions(pipeline: 'nf-core/rnaseq')
        def context = new LaunchCommandImpl.LaunchContext(
            apiEndpoint: 'https://api.cloud.seqera.io',
            userName: 'testuser',
            orgName: 'TestOrg',
            workspaceName: 'TestWS'
        )

        when:
        def result = cmd.extractLaunchResult(response, null, options, 'https://github.com/nf-core/rnaseq', context)

        then:
        result.trackingUrl == 'https://cloud.seqera.io/orgs/TestOrg/workspaces/TestWS/watch/wf-123/'
    }

    // ===== Parameter Value Parsing Tests =====

    def 'should parse parameter values correctly'() {
        given:
        def cmd = new LaunchCommandImpl()

        expect:
        cmd.parseParamValue('true') == Boolean.TRUE
        cmd.parseParamValue('false') == Boolean.FALSE
        cmd.parseParamValue('TRUE') == Boolean.TRUE
        cmd.parseParamValue('FALSE') == Boolean.FALSE
        cmd.parseParamValue('42') == 42
        cmd.parseParamValue('3.14') == 3.14
        cmd.parseParamValue('text') == 'text'
        cmd.parseParamValue(null) == null
    }

    def 'should convert kebab case to camel case'() {
        given:
        def cmd = new LaunchCommandImpl()

        expect:
        cmd.kebabToCamelCase('max-memory') == 'maxMemory'
        cmd.kebabToCamelCase('output-dir') == 'outputDir'
        cmd.kebabToCamelCase('simple') == 'simple'
        cmd.kebabToCamelCase('very-long-param-name') == 'veryLongParamName'
    }

    // ===== Utility Tests =====

    def 'should get web URL from API endpoint'() {
        given:
        def cmd = new LaunchCommandImpl()

        expect:
        cmd.getWebUrlFromApiEndpoint('https://api.cloud.seqera.io') == 'https://cloud.seqera.io'
        cmd.getWebUrlFromApiEndpoint('https://cloud.seqera.io/api') == 'https://cloud.seqera.io'
        cmd.getWebUrlFromApiEndpoint('https://custom.example.com') == 'https://custom.example.com'
    }

    // ===== Data Class Tests =====

    def 'should create LaunchContext with all fields'() {
        when:
        def context = new LaunchCommandImpl.LaunchContext(
            accessToken: 'token',
            apiEndpoint: 'https://api.cloud.seqera.io',
            userName: 'testuser',
            workspaceId: 12345L,
            orgName: 'TestOrg',
            workspaceName: 'TestWS',
            computeEnvId: 'ce-123',
            computeEnvName: 'test-ce',
            workDir: 's3://bucket/work'
        )

        then:
        context.accessToken == 'token'
        context.apiEndpoint == 'https://api.cloud.seqera.io'
        context.userName == 'testuser'
        context.workspaceId == 12345L
        context.orgName == 'TestOrg'
        context.workspaceName == 'TestWS'
        context.computeEnvId == 'ce-123'
        context.computeEnvName == 'test-ce'
        context.workDir == 's3://bucket/work'
    }

    def 'should create WorkflowLaunchResult with all fields'() {
        when:
        def result = new LaunchCommandImpl.WorkflowLaunchResult(
            workflowId: 'wf-123',
            runName: 'test-run',
            commitId: 'abc123',
            revision: 'main',
            repository: 'https://github.com/org/repo',
            trackingUrl: 'https://cloud.seqera.io/watch/wf-123'
        )

        then:
        result.workflowId == 'wf-123'
        result.runName == 'test-run'
        result.commitId == 'abc123'
        result.revision == 'main'
        result.repository == 'https://github.com/org/repo'
        result.trackingUrl == 'https://cloud.seqera.io/watch/wf-123'
    }
}
