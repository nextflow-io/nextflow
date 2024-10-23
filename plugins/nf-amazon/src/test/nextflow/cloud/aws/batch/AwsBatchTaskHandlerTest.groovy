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

import java.time.Instant

import com.amazonaws.services.batch.AWSBatch
import com.amazonaws.services.batch.model.ContainerProperties
import com.amazonaws.services.batch.model.DescribeJobDefinitionsRequest
import com.amazonaws.services.batch.model.DescribeJobDefinitionsResult
import com.amazonaws.services.batch.model.DescribeJobsRequest
import com.amazonaws.services.batch.model.DescribeJobsResult
import com.amazonaws.services.batch.model.EvaluateOnExit
import com.amazonaws.services.batch.model.JobDefinition
import com.amazonaws.services.batch.model.JobDetail
import com.amazonaws.services.batch.model.KeyValuePair
import com.amazonaws.services.batch.model.RegisterJobDefinitionRequest
import com.amazonaws.services.batch.model.RegisterJobDefinitionResult
import com.amazonaws.services.batch.model.RetryStrategy
import com.amazonaws.services.batch.model.SubmitJobRequest
import com.amazonaws.services.batch.model.SubmitJobResult
import nextflow.BuildInfo
import nextflow.Session
import nextflow.cloud.aws.config.AwsConfig
import nextflow.cloud.aws.util.S3PathFactory
import nextflow.cloud.types.CloudMachineInfo
import nextflow.cloud.types.PriceModel
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.Executor
import nextflow.fusion.FusionScriptLauncher
import nextflow.processor.Architecture
import nextflow.processor.BatchContext
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.script.BaseScript
import nextflow.script.ProcessConfig
import nextflow.util.MemoryUnit
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
        def VAR_FOO = new KeyValuePair().withName('FOO').withValue('1')
        def VAR_BAR = new KeyValuePair().withName('BAR').withValue('2')
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

        req.getJobName() == 'batch-task'
        req.getJobQueue() == 'queue1'
        req.getJobDefinition() == 'job-def:1'
        req.getContainerOverrides().getResourceRequirements().find { it.type=='VCPU'}.getValue() == '4'
        req.getContainerOverrides().getResourceRequirements().find { it.type=='MEMORY'}.getValue() == '8192'
        req.getContainerOverrides().getEnvironment() == [VAR_FOO, VAR_BAR]
        req.getContainerOverrides().getCommand() == ['bash', '-c', 'something']
        req.getRetryStrategy() == new RetryStrategy()
                .withAttempts(5)
                .withEvaluateOnExit( new EvaluateOnExit().withAction('RETRY').withOnStatusReason('Host EC2*'), new EvaluateOnExit().withOnReason('*').withAction('EXIT') )

        when:
        req = handler.newSubmitRequest(task)
        then:
        1 * handler.getSubmitCommand() >> ['bash', '-c', 'something']
        1 * handler.maxSpotAttempts() >> 0
        _ * handler.getAwsOptions() >> { new AwsOptions(awsConfig: new AwsConfig(batch: [cliPath: '/bin/aws'], region: 'eu-west-1')) }
        1 * handler.getJobQueue(task) >> 'queue1'
        1 * handler.getJobDefinition(task) >> 'job-def:1'
        1 * handler.getEnvironmentVars() >> [VAR_FOO, VAR_BAR]

        req.getJobName() == 'batch-task'
        req.getJobQueue() == 'queue1'
        req.getJobDefinition() == 'job-def:1'
        req.getContainerOverrides().getResourceRequirements().find { it.type=='VCPU'}.getValue() == '4'
        req.getContainerOverrides().getResourceRequirements().find { it.type=='MEMORY'}.getValue() == '8192'
        req.getContainerOverrides().getEnvironment() == [VAR_FOO, VAR_BAR]
        req.getContainerOverrides().getCommand() == ['bash', '-c', 'something']
        req.getRetryStrategy() == null  // <-- retry is managed by NF, hence this must be null

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

        req.getJobName() == 'batch-task'
        req.getJobQueue() == 'queue1'
        req.getJobDefinition() == 'job-def:1'
        req.getContainerOverrides().getResourceRequirements().find { it.type=='VCPU'}.getValue() == '4'
        req.getContainerOverrides().getResourceRequirements().find { it.type=='MEMORY'}.getValue() == '8192'
        req.getContainerOverrides().getCommand() == ['bash', '-c', 'something']

        when:
        def req2 = handler.newSubmitRequest(task)
        then:
        1 * handler.getSubmitCommand() >> ['bash', '-c', 'something']
        1 * handler.maxSpotAttempts() >> 5
        _ * handler.getAwsOptions() >> { new AwsOptions(awsConfig: new AwsConfig(batch: [cliPath: '/bin/aws',schedulingPriority: 9999,shareIdentifier: 'priority/high'], client:[storageEncryption: 'AES256', debug: true])) }
        1 * handler.getJobQueue(task) >> 'queue1'
        1 * handler.getJobDefinition(task) >> 'job-def:1'
        1 * handler.getEnvironmentVars() >> []

        req2.getJobName() == 'batch-task'
        req2.getJobQueue() == 'queue1'
        req2.getJobDefinition() == 'job-def:1'
        req2.getContainerOverrides().getResourceRequirements().find { it.type=='VCPU'}.getValue() == '4'
        req2.getContainerOverrides().getResourceRequirements().find { it.type=='MEMORY'}.getValue() == '8192'
        req2.getContainerOverrides().getCommand() ==['bash', '-c', 'something']
        req2.getShareIdentifier() == 'priority/high'
        req2.getSchedulingPriorityOverride() == 9999

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
        def res = req.getContainerOverrides().getResourceRequirements()
        res.size()==3
        and:
        req.getContainerOverrides().getResourceRequirements().find { it.type=='VCPU'}.getValue() == '4'
        req.getContainerOverrides().getResourceRequirements().find { it.type=='MEMORY'}.getValue() == '2048'
        req.getContainerOverrides().getResourceRequirements().find { it.type=='GPU'}.getValue() == '2'
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
        req.getJobName() == 'batch-task'
        req.getJobQueue() == 'queue1'
        req.getJobDefinition() == 'job-def:1'
        req.getTimeout() == null

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
        req.getJobName() == 'batch-task'
        req.getJobQueue() == 'queue2'
        req.getJobDefinition() == 'job-def:2'
        // minimal allowed timeout is 60 seconds
        req.getTimeout().getAttemptDurationSeconds() == 60


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
        req.getJobName() == 'batch-task'
        req.getJobQueue() == 'queue3'
        req.getJobDefinition() == 'job-def:3'
        // minimal allowed timeout is 60 seconds
        req.getTimeout().getAttemptDurationSeconds() == 3600

    }

    def 'should create an aws submit request with retry'() {

        given:
        def VAR_RETRY_MODE = new KeyValuePair().withName('AWS_RETRY_MODE').withValue('adaptive')
        def VAR_MAX_ATTEMPTS = new KeyValuePair().withName('AWS_MAX_ATTEMPTS').withValue('10')
        def VAR_METADATA_ATTEMPTS = new KeyValuePair().withName('AWS_METADATA_SERVICE_NUM_ATTEMPTS').withValue('10')
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
        req.getJobName() == 'batch-task'
        req.getJobQueue() == 'queue1'
        req.getJobDefinition() == 'job-def:1'
        // no error `retry` error strategy is defined by NF, use `maxRetries` to se Batch attempts
        req.getRetryStrategy() == new RetryStrategy()
                .withAttempts(3)
                .withEvaluateOnExit( new EvaluateOnExit().withAction('RETRY').withOnStatusReason('Host EC2*'), new EvaluateOnExit().withOnReason('*').withAction('EXIT') )
        req.getContainerOverrides().getEnvironment() == [VAR_RETRY_MODE, VAR_MAX_ATTEMPTS, VAR_METADATA_ATTEMPTS]
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
        new KeyValuePair().withName(K).withValue(V)
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

        def req = Mock(RegisterJobDefinitionRequest) {
            getJobDefinitionName() >> JOB_NAME
            getParameters() >> [ 'nf-token': JOB_ID ]
        }

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
        def client = Mock(AWSBatch)
        def handler = Spy(AwsBatchTaskHandler)
        handler.@client = client

        def req = new DescribeJobDefinitionsRequest().withJobDefinitionName(JOB_NAME)
        def res = Mock(DescribeJobDefinitionsResult)
        def job = Mock(JobDefinition)

        when:
        def result = handler.findJobDef(JOB_NAME, JOB_ID)
        then:
        1 * client.describeJobDefinitions(req) >> res
        1 * res.getJobDefinitions() >> []
        result == null

        when:
        result = handler.findJobDef(JOB_NAME, JOB_ID)
        then:
        1 * client.describeJobDefinitions(req) >> res
        1 * res.getJobDefinitions() >> [job]
        1 * job.getStatus() >> 'ACTIVE'
        1 * job.getParameters() >> ['nf-token': JOB_ID]
        1 * job.getRevision() >> 3
        result == "$JOB_NAME:3"

        when:
        result = handler.findJobDef(JOB_NAME, JOB_ID)
        then:
        1 * client.describeJobDefinitions(req) >> res
        1 * res.getJobDefinitions() >> [job]
        1 * job.getStatus() >> 'ACTIVE'
        1 * job.getParameters() >> [:]
        result == null

        when:
        result = handler.findJobDef(JOB_NAME, JOB_ID)
        then:
        1 * client.describeJobDefinitions(req) >> res
        1 * res.getJobDefinitions() >> [job]
        1 * job.getStatus() >> 'INACTIVE'
        0 * job.getParameters()
        result == null

    }

    def 'should create job definition existence' () {

        given:
        def JOB_NAME = 'foo-bar-1-0'
        def client = Mock(AWSBatch)
        def handler = Spy(AwsBatchTaskHandler)
        handler.@client = client

        def req = new RegisterJobDefinitionRequest()
        def res = Mock(RegisterJobDefinitionResult)

        when:
        def result = handler.createJobDef(req)
        then:
        1 * client.registerJobDefinition(req) >> res
        1 * res.getJobDefinitionName() >> JOB_NAME
        1 * res.getRevision() >> 10
        and:
        result == "$JOB_NAME:10"
        and:
        req.getTags().get('nextflow.io/version') == BuildInfo.version
        Instant.parse(req.getTags().get('nextflow.io/createdAt'))

    }

    def 'should add container mounts' () {
        given:
        def container = new ContainerProperties()
        def handler = Spy(AwsBatchTaskHandler)
        def mounts = [
                vol0: '/foo',
                vol1: '/foo:/bar',
                vol2: '/here:/there:ro',
                vol3: '/this:/that:rw',
        ]
        
        when:
        handler.addVolumeMountsToContainer(mounts, container)
        then:
        container.volumes.size() == 4
        container.mountPoints.size() == 4

        container.volumes[0].name == 'vol0'
        container.volumes[0].host.sourcePath == '/foo'
        container.mountPoints[0].sourceVolume == 'vol0'
        container.mountPoints[0].containerPath == '/foo'
        !container.mountPoints[0].readOnly

        container.volumes[1].name == 'vol1'
        container.volumes[1].host.sourcePath == '/foo'
        container.mountPoints[1].sourceVolume == 'vol1'
        container.mountPoints[1].containerPath == '/bar'
        !container.mountPoints[1].readOnly

        container.volumes[2].name == 'vol2'
        container.volumes[2].host.sourcePath == '/here'
        container.mountPoints[2].sourceVolume == 'vol2'
        container.mountPoints[2].containerPath == '/there'
        container.mountPoints[2].readOnly

        container.volumes[3].name == 'vol3'
        container.volumes[3].host.sourcePath == '/this'
        container.mountPoints[3].sourceVolume == 'vol3'
        container.mountPoints[3].containerPath == '/that'
        !container.mountPoints[3].readOnly
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
        result.type == 'container'
        result.parameters.'nf-token' == 'bfd3cc19ee9bdaea5b7edee94adf04bc'
        !result.containerProperties.logConfiguration
        !result.containerProperties.mountPoints
        !result.containerProperties.privileged
        
        when:
        result = handler.makeJobDefRequest(task)
        then:
        1 * handler.normalizeJobDefinitionName(IMAGE) >> JOB_NAME
        1 * handler.getAwsOptions() >> new AwsOptions(awsConfig: new AwsConfig(batch: [cliPath: '/home/conda/bin/aws', logsGroup: '/aws/batch'], region: 'us-east-1'))
        result.jobDefinitionName == JOB_NAME
        result.type == 'container'
        result.parameters.'nf-token' == 'af124f8899bcfc8a02037599f59a969a'
        result.containerProperties.logConfiguration.'LogDriver' == 'awslogs'
        result.containerProperties.logConfiguration.'Options'.'awslogs-region' == 'us-east-1'
        result.containerProperties.logConfiguration.'Options'.'awslogs-group' == '/aws/batch'
        result.containerProperties.mountPoints[0].sourceVolume == 'aws-cli'
        result.containerProperties.mountPoints[0].containerPath == '/home/conda'
        result.containerProperties.mountPoints[0].readOnly
        result.containerProperties.volumes[0].host.sourcePath == '/home/conda'
        result.containerProperties.volumes[0].name == 'aws-cli'

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
        result.type == 'container'
        result.getPlatformCapabilities() == ['FARGATE']
        result.containerProperties.getJobRoleArn() == 'the-job-role'
        result.containerProperties.getExecutionRoleArn() == 'the-exec-role'
        result.containerProperties.getResourceRequirements().find { it.type=='VCPU'}.getValue() == '1'
        result.containerProperties.getResourceRequirements().find { it.type=='MEMORY'}.getValue() == '2048'
        and:
        result.containerProperties.getEphemeralStorage().sizeInGiB == 50
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
        result.type == 'container'
        result.getPlatformCapabilities() == ['FARGATE']
        result.containerProperties.getJobRoleArn() == 'the-job-role'
        result.containerProperties.getExecutionRoleArn() == 'the-exec-role'
        result.containerProperties.getResourceRequirements().find { it.type=='VCPU'}.getValue() == '1'
        result.containerProperties.getResourceRequirements().find { it.type=='MEMORY'}.getValue() == '2048'
        and:
        result.containerProperties.getEphemeralStorage().sizeInGiB == 100
        result.containerProperties.getRuntimePlatform().getCpuArchitecture() == 'ARM64'
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
        result.type == 'container'
        result.parameters.'nf-token' == '9da434654d8c698f87da973625f57489'
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
        result.getContainerProperties().getUser() == 'foo'
        result.getContainerProperties().getPrivileged() == true

    }

    def 'should check task status' () {

        given:
        def JOB_ID = 'job-2'
        def client = Mock(AWSBatch)
        def handler = Spy(AwsBatchTaskHandler)
        handler.@client = client

        def JOB1 = new JobDetail().withJobId('job-1')
        def JOB2 = new JobDetail().withJobId('job-2')
        def JOB3 = new JobDetail().withJobId('job-3')
        def JOBS = [ JOB1, JOB2, JOB3 ]
        def resp = Mock(DescribeJobsResult)
        resp.getJobs() >> JOBS

        when:
        def result = handler.describeJob(JOB_ID)
        then:
        1 * client.describeJobs(new DescribeJobsRequest().withJobs(JOB_ID)) >> resp
        result == JOB2

    }

    def 'should check task status with empty batch collector' () {

        given:
        def collector = Mock(BatchContext)
        def JOB_ID = 'job-1'
        def client = Mock(AWSBatch)
        def handler = Spy(AwsBatchTaskHandler)
        handler.@client = client
        handler.@jobId = JOB_ID
        handler.batch(collector)

        def JOB1 = new JobDetail().withJobId('job-1')
        def JOB2 = new JobDetail().withJobId('job-2')
        def JOB3 = new JobDetail().withJobId('job-3')
        def JOBS = [ JOB1, JOB2, JOB3 ]
        def RESP = Mock(DescribeJobsResult)
        RESP.getJobs() >> JOBS

        when:
        def result = handler.describeJob(JOB_ID)
        then:
        1 * collector.contains(JOB_ID) >> false
        1 * collector.getBatchFor(JOB_ID, 100) >> ['job-1','job-2','job-3']
        1 * client.describeJobs(new DescribeJobsRequest().withJobs(['job-1','job-2','job-3'])) >> RESP
        result == JOB1

    }

    def 'should check task status with cache batch collector' () {

        given:
        def collector = Mock(BatchContext)
        def JOB_ID = 'job-1'
        def client = Mock(AWSBatch)
        def handler = Spy(AwsBatchTaskHandler)
        handler.@client = client
        handler.@jobId = JOB_ID
        handler.batch(collector)

        def JOB1 = new JobDetail().withJobId('job-1')

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
        def client = Mock(AWSBatch)
        def proxy = Mock(AwsBatchProxy)
        def handler = Spy(AwsBatchTaskHandler)
        handler.@client = proxy
        handler.task = task

        def req = Mock(SubmitJobRequest)
        def resp = Mock(SubmitJobResult)

        when:
        handler.submit()
        then:
        1 * handler.newSubmitRequest(task) >> req
        1 * handler.bypassProxy(proxy) >> client
        1 * client.submitJob(req) >> resp
        1 * resp.getJobId() >> '12345'

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
        handler.kill()
        then:
        1 * executor.shouldDeleteJob('job1') >> true
        and:
        1 * handler.terminateJob('job1') >> null

        when:
        handler.@jobId = 'job1:task2'
        handler.kill()
        then:
        1 * executor.shouldDeleteJob('job1') >> true
        and:
        1 * handler.terminateJob('job1') >> null

        when:
        handler.@jobId = 'job1:task2'
        handler.kill()
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
        def VAR_FOO = new KeyValuePair().withName('FOO').withValue('1')
        def VAR_BAR = new KeyValuePair().withName('BAR').withValue('2')
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

        req.getJobName() == 'batch-task'
        req.getJobQueue() == 'queue1'
        req.getJobDefinition() == 'job-def:1'
        req.getContainerOverrides().getResourceRequirements().find { it.type=='VCPU'}.getValue() == '4'
        req.getContainerOverrides().getResourceRequirements().find { it.type=='MEMORY'}.getValue() == '8192'
        req.getContainerOverrides().getEnvironment() == [VAR_FOO, VAR_BAR]
        req.getContainerOverrides().getCommand() == ['sh', '-c','hello']
        req.getRetryStrategy() == new RetryStrategy()
                .withAttempts(5)
                .withEvaluateOnExit( new EvaluateOnExit().withAction('RETRY').withOnStatusReason('Host EC2*'), new EvaluateOnExit().withOnReason('*').withAction('EXIT') )
        req.getTags() == [a:'b']
        req.getPropagateTags() == true
    }
    def 'get fusion submit command' () {
        given:
        def handler = Spy(AwsBatchTaskHandler) {
            fusionEnabled() >> true
            fusionLauncher() >> new FusionScriptLauncher(scheme: 's3')
            getTask() >> Mock(TaskRun) {
                getWorkDir() >> S3PathFactory.parse('s3://my-bucket/work/dir')
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
}
