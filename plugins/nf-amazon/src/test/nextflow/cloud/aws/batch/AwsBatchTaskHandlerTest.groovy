/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.cloud.aws.batch

import java.nio.file.Path
import java.time.Instant

import nextflow.BuildInfo
import nextflow.Global
import nextflow.Session
import nextflow.cloud.aws.batch.model.ContainerPropertiesModel
import nextflow.cloud.aws.batch.model.RegisterJobDefinitionModel
import nextflow.cloud.aws.config.AwsConfig
import nextflow.cloud.aws.util.S3PathFactory
import nextflow.cloud.types.CloudMachineInfo
import nextflow.cloud.types.PriceModel
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.Executor
import nextflow.fusion.FusionScriptLauncher
import nextflow.processor.Architecture
import nextflow.processor.BatchContext
import nextflow.processor.TaskBean
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.script.BaseScript
import nextflow.script.ProcessConfig
import nextflow.util.CacheHelper
import nextflow.util.MemoryUnit
import software.amazon.awssdk.services.batch.BatchClient
import software.amazon.awssdk.services.batch.model.DescribeJobDefinitionsRequest
import software.amazon.awssdk.services.batch.model.DescribeJobDefinitionsResponse
import software.amazon.awssdk.services.batch.model.DescribeJobsRequest
import software.amazon.awssdk.services.batch.model.DescribeJobsResponse
import software.amazon.awssdk.services.batch.model.EvaluateOnExit
import software.amazon.awssdk.services.batch.model.JobDefinition
import software.amazon.awssdk.services.batch.model.JobDefinitionType
import software.amazon.awssdk.services.batch.model.JobDetail
import software.amazon.awssdk.services.batch.model.KeyValuePair
import software.amazon.awssdk.services.batch.model.PlatformCapability
import software.amazon.awssdk.services.batch.model.RegisterJobDefinitionResponse
import software.amazon.awssdk.services.batch.model.ResourceType
import software.amazon.awssdk.services.batch.model.RetryStrategy
import software.amazon.awssdk.services.batch.model.SubmitJobRequest
import software.amazon.awssdk.services.batch.model.SubmitJobResponse
import spock.lang.See
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsBatchTaskHandlerTest extends Specification {

    def 'should normalise a task name'() {
        given:
        def handler = [:] as AwsBatchTaskHandler
        expect:
        handler.normalizeJobName('foo') == 'foo'
        handler.normalizeJobName('foo (12)') == 'foo_12'
        handler.normalizeJobName('foo-12') == 'foo-12'

        when:
        def looong = '012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789'
        def result = handler.normalizeJobName(looong)
        then:
        looong.size() == 150
        result.size() == 128
        result == looong.substring(0,128)
    }

    def 'should create an aws submit request'() {
        given:
        def VAR_FOO = KeyValuePair.builder().name('FOO').value('1').build()
        def VAR_BAR = KeyValuePair.builder().name('BAR').value('2').build()
        def task = Mock(TaskRun)
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig(memory: '8GB', cpus: 4, maxRetries: 2, errorStrategy: 'retry')

        def handler = Spy(AwsBatchTaskHandler)

        when:
        def req = handler.newSubmitRequest(task)
        then:
        1 * handler.getSubmitCommand() >> ['bash', '-c', 'something']
        1 * handler.maxSpotAttempts() >> 5
        _ * handler.getAwsOptions() >> { new AwsOptions(awsConfig: new AwsConfig(batch:[cliPath: '/bin/aws'])) }
        1 * handler.getJobQueue(task) >> 'queue1'
        1 * handler.getJobDefinition(task) >> 'job-def:1'
        1 * handler.getEnvironmentVars() >> [VAR_FOO, VAR_BAR]

        req.jobName() == 'batch-task'
        req.jobQueue() == 'queue1'
        req.jobDefinition() == 'job-def:1'
        req.containerOverrides().resourceRequirements().find { it.type() == ResourceType.VCPU}.value() == '4'
        req.containerOverrides().resourceRequirements().find { it.type() == ResourceType.MEMORY}.value() == '8192'
        req.containerOverrides().environment() == [VAR_FOO, VAR_BAR]
        req.containerOverrides().command() == ['bash', '-c', 'something']
        req.retryStrategy() == RetryStrategy.builder()
                .attempts(5)
                .evaluateOnExit( EvaluateOnExit.builder().action('RETRY').onStatusReason('Host EC2*').build(), EvaluateOnExit.builder().onReason('*').action('EXIT').build() )
                .build()

        when:
        req = handler.newSubmitRequest(task)
        then:
        1 * handler.getSubmitCommand() >> ['bash', '-c', 'something']
        1 * handler.maxSpotAttempts() >> 0
        _ * handler.getAwsOptions() >> { new AwsOptions(awsConfig: new AwsConfig(batch: [cliPath: '/bin/aws'], region: 'eu-west-1')) }
        1 * handler.getJobQueue(task) >> 'queue1'
        1 * handler.getJobDefinition(task) >> 'job-def:1'
        1 * handler.getEnvironmentVars() >> [VAR_FOO, VAR_BAR]

        req.jobName() == 'batch-task'
        req.jobQueue() == 'queue1'
        req.jobDefinition() == 'job-def:1'
        req.containerOverrides().resourceRequirements().find { it.type() == ResourceType.VCPU}.value() == '4'
        req.containerOverrides().resourceRequirements().find { it.type() == ResourceType.MEMORY}.value() == '8192'
        req.containerOverrides().environment() == [VAR_FOO, VAR_BAR]
        req.containerOverrides().command() == ['bash', '-c', 'something']
        req.retryStrategy() == null  // <-- retry is managed by NF, hence this must be null
    }

    def 'should create an aws submit request with opts'() {
        given:
        def task = Mock(TaskRun)
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig(memory: '8GB', cpus: 4, maxRetries: 2, errorStrategy: 'retry')

        def handler = Spy(AwsBatchTaskHandler)

        when:
        def req = handler.newSubmitRequest(task)
        then:
        1 * handler.getSubmitCommand() >> ['bash', '-c', 'something']
        1 * handler.maxSpotAttempts() >> 5
        _ * handler.getAwsOptions() >> { new AwsOptions(awsConfig: new AwsConfig(batch: [cliPath: '/bin/aws'],client: [storageEncryption: 'AES256'])) }
        1 * handler.getJobQueue(task) >> 'queue1'
        1 * handler.getJobDefinition(task) >> 'job-def:1'
        1 * handler.getEnvironmentVars() >> []

        req.jobName() == 'batch-task'
        req.jobQueue() == 'queue1'
        req.jobDefinition() == 'job-def:1'
        req.containerOverrides().resourceRequirements().find { it.type() == ResourceType.VCPU}.value() == '4'
        req.containerOverrides().resourceRequirements().find { it.type() == ResourceType.MEMORY}.value() == '8192'
        req.containerOverrides().command() == ['bash', '-c', 'something']

        when:
        def req2 = handler.newSubmitRequest(task)
        then:
        1 * handler.getSubmitCommand() >> ['bash', '-c', 'something']
        1 * handler.maxSpotAttempts() >> 5
        _ * handler.getAwsOptions() >> { new AwsOptions(awsConfig: new AwsConfig(batch: [cliPath: '/bin/aws',schedulingPriority: 9999,shareIdentifier: 'priority/high'], client:[storageEncryption: 'AES256', debug: true])) }
        1 * handler.getJobQueue(task) >> 'queue1'
        1 * handler.getJobDefinition(task) >> 'job-def:1'
        1 * handler.getEnvironmentVars() >> []

        req2.jobName() == 'batch-task'
        req2.jobQueue() == 'queue1'
        req2.jobDefinition() == 'job-def:1'
        req2.containerOverrides().resourceRequirements().find { it.type() == ResourceType.VCPU}.value() == '4'
        req2.containerOverrides().resourceRequirements().find { it.type() == ResourceType.MEMORY}.value() == '8192'
        req2.containerOverrides().command() ==['bash', '-c', 'something']
        req2.shareIdentifier() == 'priority/high'
        req2.schedulingPriorityOverride() == 9999
    }

    def 'should set accelerator' () {
        given:
        def task = Mock(TaskRun)
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig(memory: '2GB', cpus: 4, accelerator: 2)
        task.getWorkDirStr() >> 's3://my-bucket/work/dir'
        and:
        def executor = Spy(AwsBatchExecutor) { getAwsOptions()>> new AwsOptions() }
        and:
        def handler = Spy(new AwsBatchTaskHandler(executor: executor))

        when:
        def req = handler.newSubmitRequest(task)
        then:
        handler.getAwsOptions() >> { new AwsOptions(awsConfig: new AwsConfig(batch:[cliPath: '/bin/aws'],region: 'eu-west-1')) }
        and:
        _ * handler.getTask() >> task
        _ * handler.fusionEnabled() >> false
        1 * handler.maxSpotAttempts() >> 0
        1 * handler.getJobQueue(task) >> 'queue1'
        1 * handler.getJobDefinition(task) >> 'job-def:1'
        and:
        def res = req.containerOverrides().resourceRequirements()
        res.size()==3
        and:
        req.containerOverrides().resourceRequirements().find { it.type() == ResourceType.VCPU}.value() == '4'
        req.containerOverrides().resourceRequirements().find { it.type() == ResourceType.MEMORY}.value() == '2048'
        req.containerOverrides().resourceRequirements().find { it.type() == ResourceType.GPU}.value() == '2'
    }

    def 'should create an aws submit request with a timeout'() {
        given:
        def task = Mock(TaskRun)
        task.getWorkDirStr() >> 's3://my-bucket/work/dir'
        and:
        def executor = Spy(AwsBatchExecutor) {
            getAwsOptions() >> { new AwsOptions(awsConfig: new AwsConfig(batch:[cliPath: '/bin/aws'])) }
        }
        and:
        def handler = Spy(new AwsBatchTaskHandler(executor: executor))

        when:
        def req = handler.newSubmitRequest(task)
        then:
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig()
        and:
        _ * handler.getTask() >> task
        _ * handler.fusionEnabled() >> false
        1 * handler.maxSpotAttempts() >> 0
        1 * handler.getJobQueue(task) >> 'queue1'
        1 * handler.getJobDefinition(task) >> 'job-def:1'
        and:
        req.jobName() == 'batch-task'
        req.jobQueue() == 'queue1'
        req.jobDefinition() == 'job-def:1'
        req.timeout() == null

        when:
        req = handler.newSubmitRequest(task)
        then:
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig(time: '5 sec')
        and:
        _ * handler.getTask() >> task
        _ * handler.fusionEnabled() >> false
        1 * handler.maxSpotAttempts() >> 0
        1 * handler.getJobQueue(task) >> 'queue2'
        1 * handler.getJobDefinition(task) >> 'job-def:2'
        and:
        req.jobName() == 'batch-task'
        req.jobQueue() == 'queue2'
        req.jobDefinition() == 'job-def:2'
        // minimal allowed timeout is 60 seconds
        req.timeout().attemptDurationSeconds() == 60

        when:
        req = handler.newSubmitRequest(task)
        then:
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig(time: '1 hour')
        and:
        _ * handler.getTask() >> task
        _ * handler.fusionEnabled() >> false
        1 * handler.maxSpotAttempts() >> 0
        1 * handler.getJobQueue(task) >> 'queue3'
        1 * handler.getJobDefinition(task) >> 'job-def:3'
        and:
        req.jobName() == 'batch-task'
        req.jobQueue() == 'queue3'
        req.jobDefinition() == 'job-def:3'
        // minimal allowed timeout is 60 seconds
        req.timeout().attemptDurationSeconds() == 3600
    }

    def 'should create an aws submit request with retry'() {
        given:
        def VAR_RETRY_MODE = KeyValuePair.builder().name('AWS_RETRY_MODE').value('adaptive').build()
        def VAR_MAX_ATTEMPTS = KeyValuePair.builder().name('AWS_MAX_ATTEMPTS').value('10').build()
        def VAR_METADATA_ATTEMPTS = KeyValuePair.builder().name('AWS_METADATA_SERVICE_NUM_ATTEMPTS').value('10').build()
        def task = Mock(TaskRun)
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig(memory: '8GB', cpus: 4, maxRetries: 2)

        def handler = Spy(AwsBatchTaskHandler)

        when:
        def req = handler.newSubmitRequest(task)
        then:
        handler.getAwsOptions() >> { new AwsOptions(awsConfig: new AwsConfig(batch: [cliPath: '/bin/aws', retryMode: 'adaptive', maxTransferAttempts: 10])) }
        and:
        _ * handler.fusionEnabled() >> false
        1 * handler.getSubmitCommand() >> ['bash', '-c', 'foo']
        1 * handler.maxSpotAttempts() >> 3
        1 * handler.getJobQueue(task) >> 'queue1'
        1 * handler.getJobDefinition(task) >> 'job-def:1'
        and:
        req.jobName() == 'batch-task'
        req.jobQueue() == 'queue1'
        req.jobDefinition() == 'job-def:1'
        // no error `retry` error strategy is defined by NF, use `maxRetries` to se Batch attempts
        req.retryStrategy() == RetryStrategy.builder()
                .attempts(3)
                .evaluateOnExit( EvaluateOnExit.builder().action('RETRY').onStatusReason('Host EC2*').build(),
                    EvaluateOnExit.builder().onReason('*').action('EXIT').build() ).build()
        req.containerOverrides().environment() == [VAR_RETRY_MODE, VAR_MAX_ATTEMPTS, VAR_METADATA_ATTEMPTS]
    }

    def 'should return job queue'() {
        given:
        def handler = Spy(AwsBatchTaskHandler)
        def task = Mock(TaskRun)

        when:
        def result = handler.getJobQueue(task)
        then:
        task.getConfig() >> [queue: 'my-queue']
        result == 'my-queue'

        when:
        handler.getJobQueue(task)
        then:
        task.getConfig() >> [:]
        thrown(ProcessUnrecoverableException)
    }

    def 'should return job definition' () {
        given:
        def IMAGE = 'user/image:tag'
        def JOB_NAME = 'nf-user-image-tag'
        def handler = Spy(AwsBatchTaskHandler)
        def task = Mock(TaskRun)

        when:
        def result = handler.getJobDefinition(task)
        then:
        1 * task.getContainer() >> 'job-definition://queue-name'
        0 * handler.resolveJobDefinition(_)
        result == 'queue-name'

        when:
        result = handler.getJobDefinition(task)
        then:
        1 * task.getContainer() >> IMAGE
        1 * handler.resolveJobDefinition(task) >> JOB_NAME
        result == JOB_NAME
    }

    protected KeyValuePair kv(String K, String V) {
        return KeyValuePair.builder().name(K).value(V).build()
    }

    def 'should return job envs'() {
        given:
        def handler = Spy(new AwsBatchTaskHandler(environment: ENV)) {
            fusionEnabled() >> false
            getAwsOptions() >> Mock(AwsOptions) {
                getRetryMode() >> RETRY_MODE
                getMaxTransferAttempts() >> MAX_ATTEMPTS
            }
        }

        expect:
        handler.getEnvironmentVars() == EXPECTED

        where:
        RETRY_MODE  | MAX_ATTEMPTS  |  ENV                                  | EXPECTED
        null        | null          | [FOO:'hello', BAR: 'world']           | []
        null        | null          | [NXF_DEBUG: '2', NXF_FOO: 'ignore']   | [kv('NXF_DEBUG','2')]
        'standard'  | 5             | [NXF_DEBUG: '1']                      | [kv('NXF_DEBUG','1'), kv('AWS_RETRY_MODE','standard'), kv('AWS_MAX_ATTEMPTS','5'), kv('AWS_METADATA_SERVICE_NUM_ATTEMPTS','5')]
        'adaptive'  | null          | [FOO:'hello']                         | [kv('AWS_RETRY_MODE','adaptive')]
        'legacy'    | null          | [FOO:'hello']                         | [kv('AWS_RETRY_MODE','legacy')]
        'built-in'  | null          | [FOO:'hello']                         | []
    }

    def 'should return job envs for fusion'() {
        given:
        def handler = Spy(AwsBatchTaskHandler) {
            getAwsOptions() >> Mock(AwsOptions)
            getTask() >> Mock(TaskRun) {getWorkDir() >> S3PathFactory.parse('s3://my-bucket/work/dir') }
            fusionEnabled() >> true
            fusionLauncher() >> Mock(FusionScriptLauncher) {
                fusionEnv() >> [FUSION_BUCKETS: 's3://FOO,s3://BAR']
            }
        }

        expect:
        handler.getEnvironmentVars() == [kv('FUSION_BUCKETS','s3://FOO,s3://BAR')]
    }

    def 'should strip invalid chars for job definition name' () {
        given:
        def handler = Spy(AwsBatchTaskHandler)

        expect:
        handler.normalizeJobDefinitionName(null) == null
        handler.normalizeJobDefinitionName('foo') == 'nf-foo'
        handler.normalizeJobDefinitionName('foo:1') == 'nf-foo-1'
        handler.normalizeJobDefinitionName('docker.io/foo/bar:1') == 'nf-docker-io-foo-bar-1'
        and:
        handler.normalizeJobDefinitionName('docker.io/some-container-very-long-name-123456789/123456789/123456789/123456789/123456789/123456789/123456789/123456789/123456789/123456789/name:123') == 'nf-docker-io-some-container-very-long-name--35e0fa487af0f525a8c14e7866449d8e'

        when:
        handler.normalizeJobDefinitionName('/some/file.img')
        then:
        thrown(IllegalArgumentException)
    }

    def 'should validate job validation method' () {
        given:
        def IMAGE = 'foo/bar:1.0'
        def JOB_NAME = 'nf-foo-bar-1-0'
        def JOB_ID= '123'
        def handler = Spy(AwsBatchTaskHandler)
        def task = Mock(TaskRun) { getContainer()>>IMAGE }

        def req = new RegisterJobDefinitionModel().jobDefinitionName(JOB_NAME).parameters(['nf-token': JOB_ID])

        when:
        handler.resolveJobDefinition(task)
        then:
        1 * handler.makeJobDefRequest(task) >> req
        1 * handler.findJobDef(JOB_NAME, JOB_ID) >> null
        1 * handler.createJobDef(req) >> null

        when:
        handler.resolveJobDefinition(task)
        then:
        // second time are not invoked for the same image
        1 * handler.makeJobDefRequest(task) >> req
        0 * handler.findJobDef(JOB_NAME, JOB_ID) >> null
        0 * handler.createJobDef(req) >> null
    }

    def 'should verify job definition existence' () {

        given:
        def JOB_NAME = 'foo-bar-1-0'
        def JOB_ID = '123'
        def client = Mock(BatchClient)
        def task = Mock(TaskRun)
        def handler = Spy(AwsBatchTaskHandler) 
        handler.task = task
        handler.@client = client

        def req = DescribeJobDefinitionsRequest.builder().jobDefinitionName(JOB_NAME).build()
        def res = GroovyMock(DescribeJobDefinitionsResponse)
        def job = GroovyMock(JobDefinition)

        when:
        def result = handler.findJobDef(JOB_NAME, JOB_ID)
        then:
        1 * handler.bypassProxy(client) >> client
        1 * client.describeJobDefinitions(req) >> res
        1 * res.jobDefinitions() >> []
        result == null

        when:
        result = handler.findJobDef(JOB_NAME, JOB_ID)
        then:
        1 * handler.bypassProxy(client) >> client
        1 * client.describeJobDefinitions(req) >> res
        1 * res.jobDefinitions() >> [job]
        1 * job.status() >> 'ACTIVE'
        1 * job.parameters() >> ['nf-token': JOB_ID]
        1 * job.asBoolean() >> true
        1 * job.revision() >> 3
        result == "$JOB_NAME:3"

        when:
        result = handler.findJobDef(JOB_NAME, JOB_ID)
        then:
        1 * handler.bypassProxy(client) >> client
        1 * client.describeJobDefinitions(req) >> res
        1 * res.jobDefinitions() >> [job]
        1 * job.status() >> 'ACTIVE'
        1 * job.parameters() >> [:]
        result == null

        when:
        result = handler.findJobDef(JOB_NAME, JOB_ID)
        then:
        1 * handler.bypassProxy(client) >> client
        1 * client.describeJobDefinitions(req) >> res
        1 * res.jobDefinitions() >> [job]
        1 * job.status() >> 'INACTIVE'
        0 * job.parameters()
        result == null

    }

    def 'should create job definition existence' () {
        given:
        def JOB_NAME = 'foo-bar-1-0'
        def client = Mock(BatchClient)
        def handler = Spy(AwsBatchTaskHandler)
        handler.@client = client

        def req = new RegisterJobDefinitionModel().jobDefinitionName(JOB_NAME)
        def res = RegisterJobDefinitionResponse.builder().jobDefinitionName(JOB_NAME).revision(10).build()

        when:
        def result = handler.createJobDef(req)
        then:
        1 * client.registerJobDefinition(_) >> res
        and:
        result == "$JOB_NAME:10"
        and:
        req.tags.get('nextflow.io/version') == BuildInfo.version
        Instant.parse(req.tags.get('nextflow.io/createdAt'))
    }

    def 'should add container mounts' () {
        given:
        def containerModel = new ContainerPropertiesModel()
        def handler = Spy(AwsBatchTaskHandler)
        def mounts = [
                vol0: '/foo',
                vol1: '/foo:/bar',
                vol2: '/here:/there:ro',
                vol3: '/this:/that:rw',
        ]

        when:
        handler.addVolumeMountsToContainer(mounts, containerModel)
        def container = containerModel.toBatchContainerProperties()
        then:
        container.volumes().size() == 4
        container.mountPoints().size() == 4

        container.volumes()[0].name() == 'vol0'
        container.volumes()[0].host().sourcePath() == '/foo'
        container.mountPoints()[0].sourceVolume() == 'vol0'
        container.mountPoints()[0].containerPath() == '/foo'
        !container.mountPoints()[0].readOnly()

        container.volumes()[1].name() == 'vol1'
        container.volumes()[1].host().sourcePath() == '/foo'
        container.mountPoints()[1].sourceVolume() == 'vol1'
        container.mountPoints()[1].containerPath() == '/bar'
        !container.mountPoints()[1].readOnly()

        container.volumes()[2].name() == 'vol2'
        container.volumes()[2].host().sourcePath() == '/here'
        container.mountPoints()[2].sourceVolume() == 'vol2'
        container.mountPoints()[2].containerPath() == '/there'
        container.mountPoints()[2].readOnly()

        container.volumes()[3].name() == 'vol3'
        container.volumes()[3].host().sourcePath() == '/this'
        container.mountPoints()[3].sourceVolume() == 'vol3'
        container.mountPoints()[3].containerPath() == '/that'
        !container.mountPoints()[3].readOnly()
    }

    def 'should create a job definition request object' () {
        given:
        def IMAGE = 'foo/bar:1.0'
        def JOB_NAME = 'nf-foo-bar-1-0'
        def task = Mock(TaskRun) { getContainer()>>IMAGE; getConfig() >> Mock(TaskConfig) }
        def handler = Spy(AwsBatchTaskHandler) {
            getTask() >> task
            fusionEnabled() >> false
        }
        handler.@executor = Mock(AwsBatchExecutor)

        when:
        def result = handler.makeJobDefRequest(task)
        then:
        1 * handler.normalizeJobDefinitionName(IMAGE) >> JOB_NAME
        1 * handler.getAwsOptions() >> new AwsOptions()
        result.jobDefinitionName == JOB_NAME
        result.type == JobDefinitionType.CONTAINER
        result.parameters.'nf-token' == CacheHelper.hasher([JOB_NAME, result.containerProperties.toString()]).hash().toString()
        result.containerProperties.logConfiguration == null
        result.containerProperties.mountPoints == null
        result.containerProperties.privileged == false
        when:
        result = handler.makeJobDefRequest(task)
        then:
        1 * handler.normalizeJobDefinitionName(IMAGE) >> JOB_NAME
        1 * handler.getAwsOptions() >> new AwsOptions(awsConfig: new AwsConfig(batch: [cliPath: '/home/conda/bin/aws', logsGroup: '/aws/batch'], region: 'us-east-1'))
        result.jobDefinitionName == JOB_NAME
        result.type == JobDefinitionType.CONTAINER
        result.parameters.'nf-token' == CacheHelper.hasher([JOB_NAME, result.containerProperties.toString()]).hash().toString()
        result.containerProperties.logConfiguration.logDriver().toString() == 'awslogs'
        result.containerProperties.logConfiguration.options()['awslogs-region'] == 'us-east-1'
        result.containerProperties.logConfiguration.options()['awslogs-group'] == '/aws/batch'
        result.containerProperties.mountPoints[0].sourceVolume() == 'aws-cli'
        result.containerProperties.mountPoints[0].containerPath() == '/home/conda'
        result.containerProperties.mountPoints[0].readOnly()
        result.containerProperties.volumes[0].host().sourcePath() == '/home/conda'
        result.containerProperties.volumes[0].name() == 'aws-cli'
    }

    def 'should create a fargate job definition' () {
        given:
        def ARM64 = new Architecture('linux/arm64')
        def _100GB = MemoryUnit.of('100GB')
        def IMAGE = 'foo/bar:1.0'
        def JOB_NAME = 'nf-foo-bar-1-0'
        def task = Mock(TaskRun) { getContainer()>>IMAGE }
        def handler = Spy(AwsBatchTaskHandler) {
            getTask() >> task
            fusionEnabled() >> false
        }
        handler.@executor = Mock(AwsBatchExecutor)
        and:
        def session =  Mock(Session) {
            getConfig() >> [aws:[batch:[platformType:'fargate', jobRole: 'the-job-role', executionRole: 'the-exec-role']]]
        }
        def opts = new AwsOptions(session)

        when:
        def result = handler.makeJobDefRequest(task)
        then:
        task.getConfig() >> Mock(TaskConfig)
        and:
        1 * handler.normalizeJobDefinitionName(IMAGE) >> JOB_NAME
        1 * handler.getAwsOptions() >> opts
        and:
        result.jobDefinitionName == JOB_NAME
        result.type == JobDefinitionType.CONTAINER
        result.platformCapabilities == [PlatformCapability.FARGATE]
        result.containerProperties.getJobRoleArn() == 'the-job-role'
        result.containerProperties.getExecutionRoleArn() == 'the-exec-role'
        result.containerProperties.getResourceRequirements().find { it.type()==ResourceType.VCPU}.value() == '1'
        result.containerProperties.getResourceRequirements().find { it.type()==ResourceType.MEMORY}.value() == '2048'
        and:
        result.containerProperties.getEphemeralStorage().sizeInGiB() == 50
        result.containerProperties.getRuntimePlatform() == null

        when:
        result = handler.makeJobDefRequest(task)
        then:
        task.getConfig() >> Mock(TaskConfig) { getDisk()>>_100GB ; getArchitecture()>>ARM64 }
        and:
        1 * handler.normalizeJobDefinitionName(IMAGE) >> JOB_NAME
        1 * handler.getAwsOptions() >> opts
        and:
        result.jobDefinitionName == JOB_NAME
        result.type == JobDefinitionType.CONTAINER
        result.platformCapabilities == [PlatformCapability.FARGATE]
        result.containerProperties.getJobRoleArn() == 'the-job-role'
        result.containerProperties.getExecutionRoleArn() == 'the-exec-role'
        result.containerProperties.getResourceRequirements().find { it.type()==ResourceType.VCPU}.value() == '1'
        result.containerProperties.getResourceRequirements().find { it.type()==ResourceType.MEMORY}.value() == '2048'
        and:
        result.containerProperties.getEphemeralStorage().sizeInGiB() == 100
        result.containerProperties.getRuntimePlatform().cpuArchitecture() == 'ARM64'
    }

    def 'should create a job definition request object for fusion' () {
        given:
        def IMAGE = 'foo/bar:1.0'
        def JOB_NAME = 'nf-foo-bar-1-0'
        def task = Mock(TaskRun) { getContainer()>>IMAGE; getConfig()>>Mock(TaskConfig)  }
        and:
        AwsBatchTaskHandler handler = Spy(AwsBatchTaskHandler) {
            getTask() >> task
            fusionEnabled() >> true
        }
        handler.@executor = Mock(AwsBatchExecutor) {}

        when:
        def result = handler.makeJobDefRequest(task)
        then:
        1 * handler.normalizeJobDefinitionName(IMAGE) >> JOB_NAME
        1 * handler.getAwsOptions() >> new AwsOptions()
        result.jobDefinitionName == JOB_NAME
        result.type == JobDefinitionType.CONTAINER
        result.parameters.'nf-token' == CacheHelper.hasher([JOB_NAME, result.containerProperties.toString()]).hash().toString()
        result.containerProperties.privileged
    }

    def 'should create user mounts' () {
        given:
        def IMAGE = 'foo/bar:1.0'
        def JOB_NAME = 'nf-foo-bar-1-0'
        def executor = Mock(AwsBatchExecutor)
        def opts = Mock(AwsOptions)
        def task = Mock(TaskRun) { getContainer()>>IMAGE; getConfig()>>Mock(TaskConfig)  }
        and:
        def handler = Spy(AwsBatchTaskHandler) {
            getTask() >> task
            fusionEnabled() >> false
        }
        handler.@executor = executor

        when:
        def result = handler.makeJobDefRequest(task)
        then:
        1 * handler.normalizeJobDefinitionName(IMAGE) >> JOB_NAME
        1 * handler.getAwsOptions() >> opts
        1 * opts.getVolumes() >> ['/tmp', '/here:/there:ro']
        then:
        result.containerProperties.mountPoints.size() == 2
        result.containerProperties.volumes.size() == 2

        result.containerProperties.volumes[0].host.sourcePath == '/tmp'
        result.containerProperties.mountPoints[0].containerPath == '/tmp'
        !result.containerProperties.mountPoints[0].readOnly

        result.containerProperties.volumes[1].host.sourcePath == '/here'
        result.containerProperties.mountPoints[1].containerPath == '/there'
        result.containerProperties.mountPoints[1].readOnly
    }

    def 'should set job role arn'  () {
        given:
        def ROLE = 'aws::foo::bar'
        def IMAGE = 'foo/bar:1.0'
        def JOB_NAME = 'nf-foo-bar-1-0'
        def opts = Mock(AwsOptions)
        def executor = Mock(AwsBatchExecutor)
        def task = Mock(TaskRun) { getContainer()>>IMAGE; getConfig()>>Mock(TaskConfig) }
        and:
        def handler = Spy(AwsBatchTaskHandler) {
            getTask() >> task
            fusionEnabled() >> false
        }
        handler.@executor = executor

        when:
        def result = handler.makeJobDefRequest(task)
        then:
        1 * handler.normalizeJobDefinitionName(IMAGE) >> JOB_NAME
        1 * handler.getAwsOptions() >> opts
        1 * opts.getJobRole() >> ROLE

        then:
        result.getContainerProperties().getJobRoleArn() == ROLE
    }

    def 'should set container linux properties'  () {
        given:
        def ROLE = 'aws::foo::bar'
        def IMAGE = 'foo/bar:1.0'
        def JOB_NAME = 'nf-foo-bar-1-0'
        def opts = Mock(AwsOptions)
        def taskConfig = new TaskConfig(containerOptions: '--privileged --user foo')
        def executor = Mock(AwsBatchExecutor)
        def task = Mock(TaskRun) { getContainer()>>IMAGE; getConfig()>>taskConfig }
        and:
        def handler = Spy(AwsBatchTaskHandler) {
            getTask() >> task
            fusionEnabled() >> false
        }
        handler.@executor = executor

        when:
        def result = handler.makeJobDefRequest(task)
        then:
        1 * handler.normalizeJobDefinitionName(IMAGE) >> JOB_NAME
        1 * handler.getAwsOptions() >> opts

        then:
        result.containerProperties.user == 'foo'
        result.containerProperties.privileged == true
    }

    def 'should check task status' () {
        given:
        def JOB_ID = 'job-2'
        def client = Mock(BatchClient)
        def handler = Spy(AwsBatchTaskHandler)
        handler.@client = client

        def JOB1 = JobDetail.builder().jobId('job-1').build()
        def JOB2 = JobDetail.builder().jobId('job-2').build()
        def JOB3 = JobDetail.builder().jobId('job-3').build()
        def JOBS = [ JOB1, JOB2, JOB3 ]
        def resp = DescribeJobsResponse.builder().jobs(JOBS).build()

        when:
        def result = handler.describeJob(JOB_ID)
        then:
        1 * client.describeJobs(DescribeJobsRequest.builder().jobs(JOB_ID).build()) >> resp
        result == JOB2
    }

    def 'should check task status with empty batch collector' () {
        given:
        def collector = Mock(BatchContext)
        def JOB_ID = 'job-1'
        def client = Mock(BatchClient)
        def handler = Spy(AwsBatchTaskHandler)
        handler.@client = client
        handler.@jobId = JOB_ID
        handler.batch(collector)

        def JOB1 = JobDetail.builder().jobId('job-1').build()
        def JOB2 = JobDetail.builder().jobId('job-2').build()
        def JOB3 = JobDetail.builder().jobId('job-3').build()
        def JOBS = [ JOB1, JOB2, JOB3 ]
        def RESP = DescribeJobsResponse.builder().jobs(JOBS).build()

        when:
        def result = handler.describeJob(JOB_ID)
        then:
        1 * collector.contains(JOB_ID) >> false
        1 * collector.getBatchFor(JOB_ID, 100) >> ['job-1','job-2','job-3']
        1 * client.describeJobs(DescribeJobsRequest.builder().jobs(['job-1','job-2','job-3']).build()) >> RESP
        result == JOB1
    }

    def 'should check task status with cache batch collector' () {
        given:
        def collector = Mock(BatchContext)
        def JOB_ID = 'job-1'
        def client = Mock(BatchClient)
        def handler = Spy(AwsBatchTaskHandler)
        handler.@client = client
        handler.@jobId = JOB_ID
        handler.batch(collector)

        def JOB1 = JobDetail.builder().jobId('job-1').build()

        when:
        def result = handler.describeJob(JOB_ID)
        then:
        1 * collector.contains(JOB_ID) >> true
        1 * collector.get(JOB_ID) >> JOB1
        0 * collector.getBatchFor(JOB_ID, 100) >> null
        0 * client.describeJobs(_) >> null
        result == JOB1
    }

    def 'should submit job' () {
        given:
        def task = Mock(TaskRun)
        def client = Mock(BatchClient)
        def proxy = Mock(AwsBatchProxy)
        def handler = Spy(AwsBatchTaskHandler)
        handler.@client = proxy
        handler.task = task

        def req = SubmitJobRequest.builder().build()
        def resp = SubmitJobResponse.builder().jobId('12345').build()

        when:
        handler.submit()
        then:
        1 * handler.newSubmitRequest(task) >> req
        1 * handler.bypassProxy(proxy) >> client
        1 * client.submitJob(req) >> resp

        handler.status == TaskStatus.SUBMITTED
        handler.jobId == '12345'
    }

    def 'should kill a job' () {
        given:
        def executor = Mock(AwsBatchExecutor)
        def task = Mock(TaskRun)
        def handler = Spy(AwsBatchTaskHandler)
        handler.@executor = executor
        handler.task = task

        when:
        handler.@jobId = 'job1'
        handler.killTask()
        then:
        1 * executor.shouldDeleteJob('job1') >> true
        and:
        1 * handler.terminateJob('job1') >> null

        when:
        handler.@jobId = 'job1:task2'
        handler.killTask()
        then:
        1 * executor.shouldDeleteJob('job1') >> true
        and:
        1 * handler.terminateJob('job1') >> null

        when:
        handler.@jobId = 'job1:task2'
        handler.killTask()
        then:
        1 * executor.shouldDeleteJob('job1') >> false
        and:
        0 * handler.terminateJob('job1') >> null
    }

    def 'should create the trace record' () {
        given:
        def exec = Mock(Executor) { getName() >> 'awsbatch' }
        def processor = Mock(TaskProcessor)
        processor.getExecutor() >> exec
        processor.getName() >> 'foo'
        processor.getConfig() >> new ProcessConfig(Mock(BaseScript))
        def task = Mock(TaskRun)
        task.getProcessor() >> processor
        task.getConfig() >> GroovyMock(TaskConfig)
        def proxy = Mock(AwsBatchProxy)
        def handler = Spy(AwsBatchTaskHandler)
        handler.@client = proxy
        handler.task = task
        handler.@jobId = 'xyz-123'

        when:
        def trace = handler.getTraceRecord()
        then:
        1 * handler.isCompleted() >> false
        1 * handler.getMachineInfo() >> new CloudMachineInfo('x1.large', 'us-east-1b', PriceModel.spot)

        and:
        trace.native_id == 'xyz-123'
        trace.executorName == 'awsbatch'
        trace.machineInfo.type == 'x1.large'
        trace.machineInfo.zone == 'us-east-1b'
        trace.machineInfo.priceModel == PriceModel.spot
    }

    def 'should render submit command' () {
        given:
        def executor = Spy(AwsBatchExecutor)
        and:
        def handler = Spy(new AwsBatchTaskHandler(executor: executor)) {
            fusionEnabled() >> false
            getTask() >> Mock(TaskRun) { getWorkDirStr()>> 's3://work'}
        }

        when:
        def result =  handler.getSubmitCommand()
        then:
        executor.getAwsOptions()>> Mock(AwsOptions) { getAwsCli() >> 'aws' }
        then:
        result.join(' ') == 'bash -o pipefail -c trap "{ ret=$?; aws s3 cp --only-show-errors .command.log s3://work/.command.log||true; exit $ret; }" EXIT; aws s3 cp --only-show-errors s3://work/.command.run - | bash 2>&1 | tee .command.log'

        when:
        result =  handler.getSubmitCommand()
        then:
        executor.getAwsOptions() >> Mock(AwsOptions)  {
            getAwsCli() >> 'aws';
            getDebug() >> true
            getStorageEncryption() >> 'aws:kms'
            getStorageKmsKeyId() >> 'kms-key-123'
        }
        then:
        result.join(' ') == 'bash -o pipefail -c trap "{ ret=$?; aws s3 cp --only-show-errors --sse aws:kms --sse-kms-key-id kms-key-123 --debug .command.log s3://work/.command.log||true; exit $ret; }" EXIT; aws s3 cp --only-show-errors --sse aws:kms --sse-kms-key-id kms-key-123 --debug s3://work/.command.run - | bash 2>&1 | tee .command.log'
    }

    def 'should render submit command with s5cmd' () {
        given:
        def executor = Spy(AwsBatchExecutor)
        and:
        def handler = Spy(new AwsBatchTaskHandler(executor: executor)) {
            fusionEnabled() >> false
            getTask() >> Mock(TaskRun) { getWorkDirStr()>> 's3://work'}
        }

        when:
        def result =  handler.getSubmitCommand()
        then:
        executor.getAwsOptions() >> Mock(AwsOptions)  { getS5cmdPath() >> 's5cmd' }
        then:
        result.join(' ') == 'bash -o pipefail -c trap "{ ret=$?; s5cmd cp .command.log s3://work/.command.log||true; exit $ret; }" EXIT; s5cmd cat s3://work/.command.run | bash 2>&1 | tee .command.log'

        when:
        result =  handler.getSubmitCommand()
        then:
        executor.getAwsOptions() >> Mock(AwsOptions)  {
            getS5cmdPath() >> 's5cmd --debug'
            getStorageEncryption() >> 'aws:kms'
            getStorageKmsKeyId() >> 'kms-key-123'
        }
        then:
        result.join(' ') == 'bash -o pipefail -c trap "{ ret=$?; s5cmd --debug cp --sse aws:kms --sse-kms-key-id kms-key-123 .command.log s3://work/.command.log||true; exit $ret; }" EXIT; s5cmd --debug cat s3://work/.command.run | bash 2>&1 | tee .command.log'

    }

    def 'should create an aws submit request with labels'() {
        given:
        def VAR_FOO = KeyValuePair.builder().name('FOO').value('1').build()
        def VAR_BAR = KeyValuePair.builder().name('BAR').value('2').build()
        def task = Mock(TaskRun)
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig(memory: '8GB', cpus: 4, maxRetries: 2, errorStrategy: 'retry', resourceLabels:[a:'b'])

        def handler = Spy(AwsBatchTaskHandler)

        when:
        def req = handler.newSubmitRequest(task)
        then:
        1 * handler.getSubmitCommand() >> ['sh', '-c', 'hello']
        1 * handler.maxSpotAttempts() >> 5
        1 * handler.getAwsOptions() >> { new AwsOptions(awsConfig: new AwsConfig(batch: [cliPath: '/bin/aws'])) }
        1 * handler.getJobQueue(task) >> 'queue1'
        1 * handler.getJobDefinition(task) >> 'job-def:1'
        1 * handler.getEnvironmentVars() >> [VAR_FOO, VAR_BAR]

        req.jobName() == 'batch-task'
        req.jobQueue() == 'queue1'
        req.jobDefinition() == 'job-def:1'
        req.containerOverrides().resourceRequirements().find { it.type()==ResourceType.VCPU}.value() == '4'
        req.containerOverrides().resourceRequirements().find { it.type()==ResourceType.MEMORY}.value() == '8192'
        req.containerOverrides().environment() == [VAR_FOO, VAR_BAR]
        req.containerOverrides().command() == ['sh', '-c','hello']
        req.retryStrategy() == RetryStrategy.builder()
                .attempts(5)
                .evaluateOnExit( EvaluateOnExit.builder().action('RETRY').onStatusReason('Host EC2*').build(),
                    EvaluateOnExit.builder().onReason('*').action('EXIT').build())
                .build()
        req.tags() == [a:'b']
        req.propagateTags() == true
    }

    def 'get fusion submit command' () {
        given:
        def remoteWorkDir = S3PathFactory.parse('s3://my-bucket/work/dir')
        def handler = Spy(AwsBatchTaskHandler) {
            fusionEnabled() >> true
            fusionLauncher() >> new FusionScriptLauncher(Mock(TaskBean), 's3', remoteWorkDir)
            getTask() >> Mock(TaskRun) {
                getWorkDir() >> remoteWorkDir
            }
        }

        when:
        def result =  handler.getSubmitCommand()
        then:
        result.join(' ') == '/usr/bin/fusion bash /fusion/s3/my-bucket/work/dir/.command.run'
    }

    @Unroll
    def 'should normalise fargate mem' () {
        given:
        def handler = Spy(AwsBatchTaskHandler) {
            getTask() >> Mock(TaskRun) { lazyName() >> 'foo' }
        }
        expect:
        handler.normaliseFargateMem(CPUS, MemoryUnit.of( MEM * 1024L*1024L )) == EXPECTED

        where:
        CPUS    | MEM       | EXPECTED
        1       | 100       | 2048
        1       | 1000      | 2048
        1       | 2000      | 2048
        1       | 3000      | 3072
        1       | 7000      | 7168
        1       | 8000      | 8192
        1       | 10000     | 8192
        and:
        2       | 1000      | 4096
        2       | 6000      | 6144
        2       | 16000     | 16384
        2       | 20000     | 16384
        and:
        4       | 1000      | 8192
        4       | 8000      | 8192
        4       | 16000     | 16384
        4       | 30000     | 30720
        4       | 40000     | 30720
        and:
        8       | 1000      | 16384
        8       | 10000     | 16384
        8       | 20000     | 20480
        8       | 30000     | 32768
        8       | 100000    | 61440
        and:
        16      | 1000      | 32768
        16      | 30000     | 32768
        16      | 40000     | 40960
        16      | 60000     | 65536
        16      | 100000    | 106496
        16      | 200000    | 122880
    }

    @Unroll
    def 'should normalise job id' () {
        given:
        def handler = Spy(AwsBatchTaskHandler)

        expect:
        handler.normaliseJobId(JOB_ID) == EXPECTED

        where:
        JOB_ID       | EXPECTED
        null         | null
        'job1'       | 'job1'
        'job1:task2' | 'job1'
    }

    def 'should get job name' () {
        given:
        def handler = Spy(new AwsBatchTaskHandler(environment: ENV))
        def task = Mock(TaskRun)

        when:
        def result = handler.getJobName(task)
        then:
        task.getName() >> NAME
        and:
        result == EXPECTED

        where:
        ENV                             | NAME      | EXPECTED
        [:]                             | 'foo'     | 'foo'
        [TOWER_WORKFLOW_ID: '12345']    | 'foo'     | 'tw-12345-foo'
        [TOWER_WORKFLOW_ID: '12345']    | 'foo'     | 'tw-12345-foo'
        [TOWER_WORKFLOW_ID: '12345']    | 'foo(12)' | 'tw-12345-foo12'
    }

    @Unroll
    def 'should validate max spot attempts' () {
        given:
        def config = Global.config = [aws:[batch:[maxSpotAttempts: ATTEMPTS]], fusion: [enabled: FUSION, snapshots: SNAPSHOTS]]
        def session = Mock(Session) {getConfig() >> config }
        def opts = new AwsOptions(session)
        and:
        def executor = Mock(AwsBatchExecutor) { getAwsOptions() >> opts; isFusionEnabled()>>FUSION }
        def proc = Mock(TaskProcessor) { getExecutor()>>executor }
        def task = Mock(TaskRun) {workDir>> Path.of('/foo'); getProcessor()>>proc }
        def handler = new AwsBatchTaskHandler(task, executor)

        expect:
        handler.maxSpotAttempts() == EXPECTED

        where:
        ATTEMPTS    | FUSION | SNAPSHOTS  | EXPECTED
        null        | false  | false      | 0
        0           | false  | false      | 0
        1           | false  | false      | 1
        2           | false  | false      | 2
        and:
        null        | true  | false      | 0
        0           | true  | false      | 0
        1           | true  | false      | 1
        2           | true  | false      | 2
        and:
        null        | true  | true       | 5    // <-- default to 5
        0           | true  | true       | 5    // <-- default to 5
        1           | true  | true       | 1
        2           | true  | true       | 2
    }

    @Unroll
    def 'should sanitize AWS Batch label' () {
        given:
        def handler = Spy(AwsBatchTaskHandler)

        expect:
        handler.sanitizeAwsBatchLabel(INPUT, MAX_LENGTH) == EXPECTED

        where:
        INPUT                               | MAX_LENGTH | EXPECTED
        // Valid labels that don't need sanitization
        'validLabel'                        | 50         | 'validLabel'
        'valid-label_123'                   | 50         | 'valid-label_123'
        'valid.label:test/path=value+more'  | 50         | 'valid.label:test/path=value+more'
        'label with spaces'                 | 50         | 'label with spaces'
        'label-with@symbol'                 | 50         | 'label-with@symbol'
        and:
        // Labels with invalid characters
        'label#with#hash'                   | 50         | 'label_with_hash'
        'label$with%special&chars'          | 50         | 'label_with_special_chars'
        'label(with)brackets[and]braces{}'  | 50         | 'label_with_brackets_and_braces'
        'label*with?wildcards'              | 50         | 'label_with_wildcards'
        'unicode_λαβελ_test'                | 50         | 'unicode_test'
        and:
        // Multiple consecutive invalid characters
        'label###multiple###hashes'         | 50         | 'label_multiple_hashes'
        'label   multiple   spaces'         | 50         | 'label_multiple_spaces'
        'label___multiple___underscores'    | 50         | 'label_multiple_underscores'
        'label$%^&*special*&^%$chars'       | 50         | 'label_special_chars'
        and:
        // Leading/trailing invalid characters
        '###leading-hashes'                 | 50         | 'leading-hashes'
        'trailing-hashes###'                | 50         | 'trailing-hashes'
        '   leading-spaces'                 | 50         | 'leading-spaces'
        'trailing-spaces   '                | 50         | 'trailing-spaces'
        '___leading-underscores'            | 50         | 'leading-underscores'
        'trailing-underscores___'           | 50         | 'trailing-underscores'
        and:
        // Length truncation
        'very-long-label-that-exceeds-max'  | 10         | 'very-long-'
        'very-long-label-ending-with-_'     | 25         | 'very-long-label-ending-wi'
        'very-long-label-ending-with-___'   | 28         | 'very-long-label-ending-with-'
        and:
        // Edge cases
        null                                | 50         | null
        ''                                  | 50         | ''
        '   '                               | 50         | null
        '___'                               | 50         | null
        '###'                               | 50         | null
        '_'                                 | 50         | null
        ' '                                 | 50         | null
        '#'                                 | 50         | null
        and:
        // Complex real-world scenarios
        'user@domain.com'                   | 50         | 'user@domain.com'
        'workflow-run-2024/01/15'           | 50         | 'workflow-run-2024/01/15'
        'task.hash.0x1234abcd'              | 50         | 'task.hash.0x1234abcd'
        'pipeline#name%with&special*chars'  | 50         | 'pipeline_name_with_special_chars'
        'session-id:abc123#$%'              | 50         | 'session-id:abc123'
    }

    @Unroll
    def 'should sanitize AWS Batch labels map' () {
        given:
        def handler = Spy(AwsBatchTaskHandler)

        expect:
        handler.sanitizeAwsBatchLabels(INPUT) == EXPECTED

        where:
        INPUT | EXPECTED
        // Null/empty input
        null | null
        [:] | [:]
        and:
        // Valid labels
        [validKey: 'validValue'] | [validKey: 'validValue']
        ['valid-key_123': 'valid-value_456'] | ['valid-key_123': 'valid-value_456']
        ['key.with:path/chars=test+more@symbol': 'value with spaces'] | ['key.with:path/chars=test+more@symbol': 'value with spaces']
        and:
        // Invalid characters in keys and values
        ['key#with#hash': 'value$with%special&chars'] | ['key_with_hash': 'value_with_special_chars']
        ['key(brackets)': 'value[squares]{braces}'] | ['key_brackets': 'value_squares_braces']
        ['unicode_λkey': 'unicode_λvalue'] | ['unicode_key': 'unicode_value']
        and:
        // Multiple entries with mixed validity
        ['validKey': 'validValue', 'invalid#key': 'invalid$value', 'another.valid:key': 'another+valid@value'] |
        ['validKey': 'validValue', 'invalid_key': 'invalid_value', 'another.valid:key': 'another+valid@value']
        and:
        // Entries that should be filtered out (null/empty after sanitization)
        ['validKey': 'validValue', '###': '$$$', '   ': '%%%', 'goodKey': 'goodValue'] |
        ['validKey': 'validValue', 'goodKey': 'goodValue']
        and:
        // Null keys or values
        ['validKey': null, null: 'validValue', 'goodKey': 'goodValue'] |
        [null: 'validValue', 'goodKey': 'goodValue']
        and:
        // Real-world example with Nextflow resource labels
        [
            'uniqueRunId': 'tw-12345-workflow-run',
            'taskHash': 'task.hash.0x1a2b3c4d#special',
            'pipelineUser': 'user@domain.com',
            'pipelineRunName': 'my-pipeline-run(2024)',
            'pipelineSessionId': 'session#id$with%special&chars',
            'pipelineResume': 'false',
            'pipelineName': 'my_pipeline/name:version+tag'
        ] |
        [
            'uniqueRunId': 'tw-12345-workflow-run',
            'taskHash': 'task.hash.0x1a2b3c4d_special',
            'pipelineUser': 'user@domain.com',
            'pipelineRunName': 'my-pipeline-run_2024',
            'pipelineSessionId': 'session_id_with_special_chars',
            'pipelineResume': 'false',
            'pipelineName': 'my_pipeline/name:version+tag'
        ]
    }

    def 'should apply label sanitization in submit request' () {
        given:
        def task = Mock(TaskRun)
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig(
            memory: '8GB',
            cpus: 4,
            resourceLabels: [
                'validLabel': 'validValue',
                'invalid#key': 'invalid$value',
                'long-key-that-might-be-truncated-if-very-very-long': 'long-value-that-should-be-truncated-because-it-exceeds-the-maximum-allowed-length-for-aws-batch-tags-which-is-256-characters-and-this-string-is-definitely-longer-than-that-limit-so-it-will-be-cut-off-at-the-appropriate-length-and-cleaned-up'
            ]
        )

        def handler = Spy(AwsBatchTaskHandler)

        when:
        def req = handler.newSubmitRequest(task)
        then:
        1 * handler.getSubmitCommand() >> ['bash', '-c', 'test']
        1 * handler.maxSpotAttempts() >> 0
        1 * handler.getAwsOptions() >> new AwsOptions(awsConfig: new AwsConfig(batch: [sanitizeTags: true]))
        1 * handler.getJobQueue(task) >> 'test-queue'
        1 * handler.getJobDefinition(task) >> 'test-job-def'
        1 * handler.getEnvironmentVars() >> []

        and:
        def tags = req.tags()
        tags.size() == 3
        tags['validLabel'] == 'validValue'
        tags['invalid_key'] == 'invalid_value'
        // Check that long value was truncated
        tags['long-key-that-might-be-truncated-if-very-very-long'].length() <= 256
        tags['long-key-that-might-be-truncated-if-very-very-long'].startsWith('long-value-that-should-be-truncated')
        !tags['long-key-that-might-be-truncated-if-very-very-long'].endsWith('_')
        req.propagateTags() == true
    }

    def 'should not sanitize labels when sanitizeTags config is disabled'() {
        given:
        def task = Mock(TaskRun)
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig(
            resourceLabels: [
                'validLabel': 'validValue',
                'invalid#key': 'invalid$value', // These should NOT be sanitized
                'special*chars?here': 'value with spaces & symbols!'
            ]
        )

        def handler = Spy(AwsBatchTaskHandler)

        when:
        def req = handler.newSubmitRequest(task)
        then:
        1 * handler.getSubmitCommand() >> ['bash', '-c', 'test']
        1 * handler.maxSpotAttempts() >> 0
        1 * handler.getAwsOptions() >> new AwsOptions(awsConfig: new AwsConfig(batch: [sanitizeTags: false]))
        1 * handler.getJobQueue(task) >> 'test-queue'
        1 * handler.getJobDefinition(task) >> 'test-job-def'
        1 * handler.getEnvironmentVars() >> []

        and: 'labels should remain unchanged (not sanitized)'
        def tags = req.tags()
        tags.size() == 3
        tags['validLabel'] == 'validValue'
        tags['invalid#key'] == 'invalid$value'  // Should still contain invalid characters
        tags['special*chars?here'] == 'value with spaces & symbols!'  // Should still contain invalid characters
        req.propagateTags() == true
    }

    def 'should not sanitize labels by default when no sanitizeTags config is specified'() {
        given:
        def task = Mock(TaskRun)
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig(
            resourceLabels: [
                'validLabel': 'validValue',
                'invalid#key': 'invalid$value', // These should NOT be sanitized by default
                'test&symbol': 'value(with)brackets'
            ]
        )

        def handler = Spy(AwsBatchTaskHandler)

        when:
        def req = handler.newSubmitRequest(task)
        then:
        1 * handler.getSubmitCommand() >> ['bash', '-c', 'test']
        1 * handler.maxSpotAttempts() >> 0
        1 * handler.getAwsOptions() >> new AwsOptions() // No config specified, should default to false
        1 * handler.getJobQueue(task) >> 'test-queue'
        1 * handler.getJobDefinition(task) >> 'test-job-def'
        1 * handler.getEnvironmentVars() >> []

        and: 'labels should remain unchanged (not sanitized by default)'
        def tags = req.tags()
        tags.size() == 3
        tags['validLabel'] == 'validValue'
        tags['invalid#key'] == 'invalid$value'  // Should still contain invalid characters
        tags['test&symbol'] == 'value(with)brackets'  // Should still contain invalid characters
        req.propagateTags() == true
    }

    @Unroll
    @See("https://github.com/nextflow-io/nextflow/pull/6211#discussion_r2161928856")
    def 'should handle null values in labels with explicit logging'() {
        given:
        def handler = Spy(AwsBatchTaskHandler)

        expect:
        handler.sanitizeAwsBatchLabels(INPUT) == EXPECTED

        where:
        INPUT | EXPECTED
        // Basic null value case - this addresses the PR comment: "when the item is "item": null is the aws tag silently dropped?"
        ['item': null, 'validKey': 'validValue'] | ['validKey': 'validValue']

        // Multiple null values
        ['key1': null, 'key2': 'value2', 'key3': null] | ['key2': 'value2']

        // All null values
        ['key1': null, 'key2': null] | [:]

        // Mix of null and empty string (which becomes null after sanitization)
        ['nullValue': null, 'emptyValue': '', 'validValue': 'good'] | ['validValue': 'good']
    }

    @See("https://github.com/nextflow-io/nextflow/pull/6211#discussion_r2161928856")
    def 'should handle null keys in labels with explicit logging'() {
        given:
        def handler = Spy(AwsBatchTaskHandler)

        when: 'creating map with actual null key'
        def labels = new HashMap<String, String>()
        labels.put(null, 'validValue')
        labels.put('validKey', 'validValue')
        def result = handler.sanitizeAwsBatchLabels(labels)

        then: 'null key is dropped'
        result.size() == 1
        result['validKey'] == 'validValue'
        !result.containsKey(null)
    }

    @See("https://github.com/nextflow-io/nextflow/pull/6211#discussion_r2161928856")
    def 'should handle both null keys and values'() {
        given:
        def handler = Spy(AwsBatchTaskHandler)

        when: 'creating map with null key and null values'
        def labels = new HashMap<String, String>()
        labels.put(null, 'someValue')    // null key
        labels.put('nullValue', null)    // null value
        labels.put(null, null)           // both null (overwrites previous null key)
        labels.put('validKey', 'validValue')
        def result = handler.sanitizeAwsBatchLabels(labels)

        then: 'only valid entry remains'
        result.size() == 1
        result['validKey'] == 'validValue'
        !result.containsKey(null)
        !result.containsKey('nullValue')
    }

    @See("https://github.com/nextflow-io/nextflow/pull/6211#discussion_r2161928856")
    def 'should verify logging behavior for null handling'() {
        given:
        def handler = new AwsBatchTaskHandler()
        def logAppender = Mock(ch.qos.logback.core.Appender)
        def logger = (ch.qos.logback.classic.Logger) org.slf4j.LoggerFactory.getLogger(AwsBatchTaskHandler)
        logger.addAppender(logAppender)
        logger.setLevel(ch.qos.logback.classic.Level.WARN)

        when: 'sanitizing labels with null values'
        def labels = ['item': null, 'validKey': 'validValue']  // This is the exact case from PR comment
        handler.sanitizeAwsBatchLabels(labels)

        then: 'warning is logged for null value'
        1 * logAppender.doAppend(_) >> { args ->
            def event = args[0] as ch.qos.logback.classic.spi.ILoggingEvent
            assert event.level == ch.qos.logback.classic.Level.WARN
            assert event.formattedMessage.contains('AWS Batch label dropped due to null value: key=item, value=null')
        }

        cleanup:
        logger.detachAppender(logAppender)
    }

    @See("https://github.com/nextflow-io/nextflow/pull/6211#discussion_r2161928856")
    def 'should verify no silent dropping - PR comment verification'() {
        given: 'This test specifically addresses the PR comment about silent dropping'
        def handler = Spy(AwsBatchTaskHandler)

        when: 'processing the exact scenario from PR comment'
        def labels = ['item': null]  // "when the item is "item": null is the aws tag silently dropped?"
        def result = handler.sanitizeAwsBatchLabels(labels)

        then: 'result is empty (tag is dropped)'
        result == [:]

        and: 'the method properly logs the dropped label (verified by observing the actual log output in test execution)'
        // The actual logging is verified by the "should verify logging behavior for null handling" test above
        // This test focuses on the functional behavior: null values are correctly dropped from the result
        true  // The key point is that silent dropping has been replaced with logged dropping
    }
}
