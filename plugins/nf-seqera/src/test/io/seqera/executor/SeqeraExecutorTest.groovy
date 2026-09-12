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

package io.seqera.executor

import io.seqera.config.ExecutorOpts
import io.seqera.config.SeqeraConfig
import io.seqera.sched.api.schema.v1a1.CreateRunRequest
import io.seqera.sched.api.schema.v1a1.CreateRunResponse
import io.seqera.sched.client.SchedClient
import io.seqera.sched.client.SchedClientConfig
import nextflow.Session
import nextflow.SysEnv
import nextflow.platform.PlatformHelper
import nextflow.script.WorkflowMetadata
import spock.lang.Specification

/**
 * Tests for SeqeraExecutor client configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SeqeraExecutorTest extends Specification {

    def cleanup() {
        SysEnv.pop()
    }

    def 'should create client config with config settings'() {
        given:
        SysEnv.push([:])

        when:
        def config = buildClientConfig(
            [endpoint: 'https://sched.example.com', region: 'us-west-2'],
            [endpoint: 'https://api.platform.example.com', accessToken: 'config-access-token', refreshToken: 'config-refresh-token']
        )

        then:
        config.endpoint == 'https://sched.example.com'
        config.platformUrl == 'https://api.platform.example.com'
        config.accessToken == 'config-access-token'
        config.refreshToken == 'config-refresh-token'
    }

    def 'should bound each attempt with the configured request timeout'() {
        given:
        SysEnv.push([:])

        when:
        def config = buildClientConfig(
            [endpoint: 'https://sched.example.com', httpClient: [requestTimeout: '5 sec']],
            [endpoint: 'https://api.platform.example.com', accessToken: 'tok']
        )

        then:
        config.requestTimeout == java.time.Duration.ofSeconds(5)
    }

    def 'should leave responses unbounded when the request timeout is zero'() {
        given:
        SysEnv.push([:])

        when:
        def config = buildClientConfig(
            [endpoint: 'https://sched.example.com', httpClient: [requestTimeout: '0 sec']],
            [endpoint: 'https://api.platform.example.com', accessToken: 'tok']
        )

        then:
        config.requestTimeout == null
    }

    def 'should bound the connect phase with the configured connect timeout'() {
        given:
        SysEnv.push([:])

        when:
        def config = buildClientConfig(
            [endpoint: 'https://sched.example.com', httpClient: [connectTimeout: '3 sec']],
            [endpoint: 'https://api.platform.example.com', accessToken: 'tok']
        )

        then:
        config.connectTimeout == java.time.Duration.ofSeconds(3)
    }

    def 'should carry the default timeouts when the http client scope is absent'() {
        given:
        SysEnv.push([:])

        when:
        def config = buildClientConfig(
            [endpoint: 'https://sched.example.com'],
            [endpoint: 'https://api.platform.example.com', accessToken: 'tok']
        )

        then:
        config.connectTimeout == java.time.Duration.ofSeconds(10)
        config.requestTimeout == java.time.Duration.ofSeconds(45)
    }

    def 'should create client config with env variable settings'() {
        given:
        SysEnv.push([
            TOWER_API_ENDPOINT: 'https://api.env.example.com',
            TOWER_ACCESS_TOKEN: 'env-access-token',
            TOWER_REFRESH_TOKEN: 'env-refresh-token'
        ])

        when:
        def config = buildClientConfig(
            [endpoint: 'https://sched.example.com', region: 'us-west-2'],
            [:]
        )

        then:
        config.endpoint == 'https://sched.example.com'
        config.platformUrl == 'https://api.env.example.com'
        config.accessToken == 'env-access-token'
        config.refreshToken == 'env-refresh-token'
    }

    def 'should use default platform url when not configured'() {
        given:
        SysEnv.push([:])

        when:
        def config = buildClientConfig(
            [endpoint: 'https://sched.example.com', region: 'us-west-2'],
            [accessToken: 'my-token']
        )

        then:
        config.endpoint == 'https://sched.example.com'
        config.platformUrl == 'https://api.cloud.seqera.io'
        config.accessToken == 'my-token'
        config.refreshToken == null
    }

    def 'should prefer config over env variables'() {
        given:
        SysEnv.push([
            TOWER_API_ENDPOINT: 'https://api.env.example.com',
            TOWER_ACCESS_TOKEN: 'env-access-token',
            TOWER_REFRESH_TOKEN: 'env-refresh-token'
        ])

        when:
        def config = buildClientConfig(
            [endpoint: 'https://sched.example.com', region: 'us-west-2'],
            [endpoint: 'https://api.config.example.com', accessToken: 'config-access-token', refreshToken: 'config-refresh-token']
        )

        then:
        config.endpoint == 'https://sched.example.com'
        config.platformUrl == 'https://api.config.example.com'
        config.accessToken == 'config-access-token'
        config.refreshToken == 'config-refresh-token'
    }

    def 'should set fusion default version when not configured' () {
        given:
        SysEnv.push([:])
        def fusionConfig = [enabled: true]
        def config = [fusion: fusionConfig]
        def session = Mock(Session) { getConfig() >> config }
        def executor = new SeqeraExecutor(session: session)

        when:
        executor.applyFusionDefaults()

        then:
        fusionConfig.targetVersion == '2.6'
    }

    def 'should not override fusion version when containerConfigUrl is set' () {
        given:
        SysEnv.push([:])
        def fusionConfig = [enabled: true, containerConfigUrl: 'https://custom.url/v3.0-amd64.json']
        def config = [fusion: fusionConfig]
        def session = Mock(Session) { getConfig() >> config }
        def executor = new SeqeraExecutor(session: session)

        when:
        executor.applyFusionDefaults()

        then:
        fusionConfig.targetVersion == null
    }

    def 'should be secret-native so Nextflow suppresses the local-store secrets snippet'() {
        given:
        SysEnv.push([:])

        expect:
        new SeqeraExecutor().isSecretNative()
    }

    def 'should build the run baseline from the auto labels and the config-level process.resourceLabels'() {
        given:
        SysEnv.push([:])
        def executor = new SeqeraExecutor()
        executor.session = Mock(Session) {
            getAutoResourceLabels() >> ['nextflow.io/runName': 'crazy_darwin']
            getConfig() >> [process: [resourceLabels: [team: 'a', priority: 7]]]
        }

        when:
        executor.computeRunResourceLabels()

        then: 'the baseline is the complete label set, so the per-task delta can collapse to empty'
        executor.runResourceLabels == ['nextflow.io/runName': 'crazy_darwin', team: 'a', priority: '7']
    }

    def 'should let the config-level process.resourceLabels win over an auto label'() {
        given:
        SysEnv.push([:])
        def executor = new SeqeraExecutor()
        executor.session = Mock(Session) {
            getAutoResourceLabels() >> ['nextflow.io/runName': 'crazy_darwin']
            getConfig() >> [process: [resourceLabels: ['nextflow.io/runName': 'custom']]]
        }

        when:
        executor.computeRunResourceLabels()

        then:
        executor.runResourceLabels == ['nextflow.io/runName': 'custom']
    }

    def 'should yield an empty run baseline when neither auto nor config-level labels are given'() {
        given:
        SysEnv.push([:])
        def executor = new SeqeraExecutor()
        executor.session = Mock(Session) {
            getAutoResourceLabels() >> [:]
            getConfig() >> [:]
        }

        when:
        executor.computeRunResourceLabels()

        then:
        executor.runResourceLabels == [:]
    }

    def 'should skip a dynamic process.resourceLabels while keeping the auto labels in the baseline'() {
        given:
        SysEnv.push([:])
        def executor = new SeqeraExecutor()
        def dynamic = { [team: 'a', priority: 7] }
        executor.session = Mock(Session) {
            getAutoResourceLabels() >> ['nextflow.io/runName': 'crazy_darwin']
            getConfig() >> [process: [resourceLabels: dynamic]]
        }

        when:
        executor.computeRunResourceLabels()

        then:
        noExceptionThrown()
        executor.runResourceLabels == ['nextflow.io/runName': 'crazy_darwin']
    }

    def 'createRun populates CreateRunRequest.labels with config-level resourceLabels merged with auto-labels'() {
        given:
        SysEnv.push([:])
        CreateRunRequest captured = null
        def mockClient = Mock(SchedClient) {
            createRun(_) >> { args ->
                captured = args[0] as CreateRunRequest
                new CreateRunResponse().runId('run-1')
            }
        }
        def autoLabels = [
            'nextflow.io/projectName': 'my-project',
            'nextflow.io/runName': 'test-run',
            'seqera.io/platform/workspaceId': '1234',
            'seqera.io/platform/computeEnvId': 'ce-abc'
        ]
        def sessionConfig = [
            process: [resourceLabels: [team: 'platform', priority: 3]],
            seqera: [executor: [endpoint: 'https://sched.example.com', provider: 'aws', region: 'us-east-1']],
            tower: [:]
        ]
        def session = Mock(Session) {
            getAutoResourceLabels() >> autoLabels
            getConfig() >> sessionConfig
            getWorkflowMetadata() >> Mock(WorkflowMetadata) { getPlatform() >> null }
            getWorkDir() >> java.nio.file.Paths.get('/work')
            getRunName() >> 'test-run'
        }
        def seqeraOpts = new ExecutorOpts(endpoint: 'https://sched.example.com', provider: 'aws', region: 'us-east-1')
        def executor = new SeqeraExecutor()
        executor.session = session
        executor.@seqeraConfig = seqeraOpts
        executor.@client = mockClient

        when:
        executor.createRun()

        then:
        captured != null
        captured.getLabels()['team'] == 'platform'
        captured.getLabels()['priority'] == '3'
        captured.getLabels()['nextflow.io/projectName'] == 'my-project'
        captured.getLabels()['nextflow.io/runName'] == 'test-run'
        captured.getLabels()['seqera.io/platform/workspaceId'] == '1234'
        captured.getLabels()['seqera.io/platform/computeEnvId'] == 'ce-abc'

        cleanup:
        executor.batchSubmitter?.shutdown()
    }

    def 'createRun honours the deprecated seqera.executor.autoLabels over tower.autoLabels'() {
        given:
        SysEnv.push([:])
        CreateRunRequest captured = null
        def mockClient = Mock(SchedClient) {
            createRun(_) >> { args ->
                captured = args[0] as CreateRunRequest
                new CreateRunResponse().runId('run-1')
            }
        }
        // a real session, so that the label set the executor sends is the one resolved by
        // `Session#getAutoResourceLabels` out of the two config options
        def session = new Session([
            seqera: [executor: [endpoint: 'https://sched.example.com', autoLabels: 'projectName']],
            tower: [autoLabels: 'runName']
        ])
        session.@workflowMetadata = Mock(WorkflowMetadata) {
            getProjectName() >> 'my-project'
            getRunName() >> 'test-run'
        }
        def executor = new SeqeraExecutor()
        executor.session = session
        executor.@seqeraConfig = new ExecutorOpts(endpoint: 'https://sched.example.com', autoLabels: 'projectName')
        executor.@client = mockClient

        when:
        executor.createRun()

        then:
        captured != null
        captured.getLabels() == ['nextflow.io/projectName': 'my-project']

        cleanup:
        executor.batchSubmitter?.shutdown()
    }

    def 'createRun passes provider, strategy and region to CreateRunRequest'() {
        given:
        SysEnv.push([:])
        CreateRunRequest captured = null
        def mockClient = Mock(SchedClient) {
            createRun(_) >> { args ->
                captured = args[0] as CreateRunRequest
                new CreateRunResponse().runId('run-1')
            }
        }
        def workflowMeta = Mock(WorkflowMetadata) {
            getPlatform() >> null
        }
        def session = Mock(Session) {
            getConfig() >> [tower: [:]]
            getWorkflowMetadata() >> workflowMeta
            getWorkDir() >> java.nio.file.Paths.get('/work')
            getRunName() >> 'test-run'
        }
        def seqeraOpts = new ExecutorOpts(
            endpoint: 'https://sched.example.com',
            provider: 'aws',
            strategy: 'vm',
            region: 'eu-west-1'
        )
        def executor = new SeqeraExecutor()
        executor.session = session
        executor.@seqeraConfig = seqeraOpts
        executor.@client = mockClient

        when:
        executor.createRun()

        then:
        captured != null
        captured.getProvider() == 'aws'
        captured.getStrategy() == 'vm'
        captured.getRegion() == 'eu-west-1'

        cleanup:
        executor.batchSubmitter?.shutdown()
    }

    def 'createRun passes maxCpusPerUser to CreateRunRequest.schedulingRequirement'() {
        given:
        SysEnv.push([:])
        CreateRunRequest captured = null
        def mockClient = Mock(SchedClient) {
            createRun(_) >> { args ->
                captured = args[0] as CreateRunRequest
                new CreateRunResponse().runId('run-1')
            }
        }
        def workflowMeta = Mock(WorkflowMetadata) {
            getPlatform() >> null
        }
        def session = Mock(Session) {
            getConfig() >> [tower: [:]]
            getWorkflowMetadata() >> workflowMeta
            getWorkDir() >> java.nio.file.Paths.get('/work')
            getRunName() >> 'test-run'
        }
        def seqeraOpts = new ExecutorOpts(
            endpoint: 'https://sched.example.com',
            schedulingRequirement: [maxCpusPerUser: 16]
        )
        def executor = new SeqeraExecutor()
        executor.session = session
        executor.@seqeraConfig = seqeraOpts
        executor.@client = mockClient

        when:
        executor.createRun()

        then:
        captured != null
        captured.getSchedulingRequirement() != null
        captured.getSchedulingRequirement().getMaxCpusPerUser() == 16

        cleanup:
        executor.batchSubmitter?.shutdown()
    }

    def 'createRun omits schedulingRequirement when maxCpusPerUser is not set'() {
        given:
        SysEnv.push([:])
        CreateRunRequest captured = null
        def mockClient = Mock(SchedClient) {
            createRun(_) >> { args ->
                captured = args[0] as CreateRunRequest
                new CreateRunResponse().runId('run-1')
            }
        }
        def workflowMeta = Mock(WorkflowMetadata) {
            getPlatform() >> null
        }
        def session = Mock(Session) {
            getConfig() >> [tower: [:]]
            getWorkflowMetadata() >> workflowMeta
            getWorkDir() >> java.nio.file.Paths.get('/work')
            getRunName() >> 'test-run'
        }
        def seqeraOpts = new ExecutorOpts(
            endpoint: 'https://sched.example.com'
        )
        def executor = new SeqeraExecutor()
        executor.session = session
        executor.@seqeraConfig = seqeraOpts
        executor.@client = mockClient

        when:
        executor.createRun()

        then:
        captured != null
        captured.getSchedulingRequirement() == null

        cleanup:
        executor.batchSubmitter?.shutdown()
    }

    def 'createRun passes providerConfig to CreateRunRequest'() {
        given:
        SysEnv.push([:])
        CreateRunRequest captured = null
        def mockClient = Mock(SchedClient) {
            createRun(_) >> { args ->
                captured = args[0] as CreateRunRequest
                new CreateRunResponse().runId('run-1')
            }
        }
        def workflowMeta = Mock(WorkflowMetadata) {
            getPlatform() >> null
        }
        def session = Mock(Session) {
            getConfig() >> [tower: [:]]
            getWorkflowMetadata() >> workflowMeta
            getWorkDir() >> java.nio.file.Paths.get('/work')
            getRunName() >> 'test-run'
        }
        def seqeraOpts = new ExecutorOpts(
            endpoint: 'https://sched.example.com',
            providerConfig: [subnetId: 'subnet-1', securityGroup: 'sg-2']
        )
        def executor = new SeqeraExecutor()
        executor.session = session
        executor.@seqeraConfig = seqeraOpts
        executor.@client = mockClient

        when:
        executor.createRun()

        then:
        captured != null
        captured.getProviderConfig() == [subnetId: 'subnet-1', securityGroup: 'sg-2']

        cleanup:
        executor.batchSubmitter?.shutdown()
    }

    def 'createRun publishes the run id to the platform metadata'() {
        given:
        SysEnv.push([:])
        def mockClient = Mock(SchedClient) {
            createRun(_) >> new CreateRunResponse().runId('run-xyz')
        }
        def platformMeta = Mock(nextflow.script.PlatformMetadata)
        def workflowMeta = Mock(WorkflowMetadata) {
            getPlatform() >> platformMeta
        }
        def session = Mock(Session) {
            getConfig() >> [tower: [:]]
            getWorkflowMetadata() >> workflowMeta
            getWorkDir() >> java.nio.file.Paths.get('/work')
            getRunName() >> 'test-run'
        }
        def seqeraOpts = new ExecutorOpts(endpoint: 'https://sched.example.com', provider: 'aws', region: 'eu-west-1')
        def executor = new SeqeraExecutor()
        executor.session = session
        executor.@seqeraConfig = seqeraOpts
        executor.@client = mockClient

        when:
        executor.createRun()

        then:
        executor.runId == 'run-xyz'
        and:
        1 * platformMeta.setSchedRunId('run-xyz')

        cleanup:
        executor.batchSubmitter?.shutdown()
    }


    /**
     * Builds a SchedClientConfig through the same code {@link SeqeraExecutor#createClient()}
     * uses, so removing a setting from that chain fails these specs rather than passing.
     */
    private SchedClientConfig buildClientConfig(Map executorOpts, Map towerConfig) {
        def seqeraConfig = new SeqeraConfig([executor: executorOpts]).executor
        return SeqeraExecutor.clientConfig(seqeraConfig, towerConfig)
    }
}
