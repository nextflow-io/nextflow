package nextflow.executor
import java.nio.file.Paths

import com.amazonaws.services.batch.AWSBatchClient
import com.amazonaws.services.batch.model.DescribeJobDefinitionsRequest
import com.amazonaws.services.batch.model.DescribeJobDefinitionsResult
import com.amazonaws.services.batch.model.JobDefinition
import com.amazonaws.services.batch.model.KeyValuePair
import com.amazonaws.services.batch.model.RegisterJobDefinitionRequest
import com.amazonaws.services.batch.model.RegisterJobDefinitionResult
import com.amazonaws.services.batch.model.RetryStrategy
import nextflow.Session
import nextflow.exception.ProcessUnrecoverableException
import nextflow.processor.TaskBean
import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
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
    }

    @Unroll
    def 'should return aws options'() {
        given:
        def cfg = [
                aws: [client: [
                        uploadStorageClass: awsStorClass,
                        storageEncryption  : awsStorEncrypt]],
                executor: [
                        awscli: awscliPath
                ]
        ]
        def session = new Session(cfg)
        def executor = Mock(AwsBatchExecutor)
        executor.getSession() >> session

        def handler = new AwsBatchTaskHandler(executor: executor)

        when:
        def opts = handler.getAwsOptions()
        then:
        opts.cliPath == awscliPath
        opts.storageClass == awsStorClass
        opts.storageEncryption == awsStorEncrypt

        where:
        awscliPath      | awsStorClass | awsStorEncrypt
        null            | null         | null
        '/foo/bin/aws'  | 'STANDARD'   | 'AES256'

    }

    def 'should validate aws options' () {

        when:
        def opts = new AwsOptions()
        then:
        opts.getCliPath() == null
        opts.getStorageClass() == null
        opts.getStorageEncryption() == null

        when:
        opts = new AwsOptions(cliPath: '/foo/bin/aws', storageClass: 'STANDARD', storageEncryption: 'AES256')
        then:
        opts.getCliPath() == '/foo/bin/aws'
        opts.getStorageClass() == 'STANDARD'
        opts.getStorageEncryption() == 'AES256'

        when:
        opts = new AwsOptions(storageClass: 'foo')
        then:
        opts.getStorageClass() == null

        when:
        opts = new AwsOptions(storageEncryption: 'abr')
        then:
        opts.getStorageEncryption() == null

        when:
        new AwsOptions(cliPath: 'bin/aws')
        then:
        thrown(ProcessUnrecoverableException)

        when:
        new AwsOptions(cliPath: '/foo/aws')
        then:
        thrown(ProcessUnrecoverableException)
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
        1 * handler.getAwsOptions() >> { new AwsOptions(cliPath: '/bin/aws') }
        1 * handler.getJobQueue(task) >> 'queue1'
        1 * handler.getJobDefinition(task) >> 'job-def:1'
        1 * handler.getEnvironmentVars() >> [VAR_FOO, VAR_BAR]
        1 * handler.wrapperFile >> Paths.get('/bucket/test/.command.run')
        1 * handler.getLogFile() >> Paths.get('/bucket/test/.command.log')

        req.getJobName() == 'batchtask'
        req.getJobQueue() == 'queue1'
        req.getJobDefinition() == 'job-def:1'
        req.getContainerOverrides().getVcpus() == 4
        req.getContainerOverrides().getMemory() == 8192
        req.getContainerOverrides().getEnvironment() == [VAR_FOO, VAR_BAR]
        req.getContainerOverrides().getCommand() == ['bash', '-o','pipefail','-c', "/bin/aws s3 cp s3://bucket/test/.command.run - | bash 2>&1 | /bin/aws s3 cp - s3://bucket/test/.command.log".toString()]
        req.getRetryStrategy() == null  // <-- retry is managed by NF, hence this must be null

        when:
        req = handler.newSubmitRequest(task)
        then:
        1 * handler.getAwsOptions() >> { new AwsOptions(cliPath: '/bin/aws', region: 'eu-west-1') }
        1 * handler.getJobQueue(task) >> 'queue1'
        1 * handler.getJobDefinition(task) >> 'job-def:1'
        1 * handler.getEnvironmentVars() >> [VAR_FOO, VAR_BAR]
        1 * handler.wrapperFile >> Paths.get('/bucket/test/.command.run')
        1 * handler.getLogFile() >> Paths.get('/bucket/test/.command.log')

        req.getJobName() == 'batchtask'
        req.getJobQueue() == 'queue1'
        req.getJobDefinition() == 'job-def:1'
        req.getContainerOverrides().getVcpus() == 4
        req.getContainerOverrides().getMemory() == 8192
        req.getContainerOverrides().getEnvironment() == [VAR_FOO, VAR_BAR]
        req.getContainerOverrides().getCommand() == ['bash', '-o','pipefail','-c', "export AWS_DEFAULT_REGION='eu-west-1'; /bin/aws s3 cp s3://bucket/test/.command.run - | bash 2>&1 | /bin/aws s3 cp - s3://bucket/test/.command.log".toString()]
        req.getRetryStrategy() == null  // <-- retry is managed by NF, hence this must be null

    }

    def 'should create an aws submit request with retry'() {

        given:
        def VAR_FOO = new KeyValuePair().withName('FOO').withValue('1')
        def VAR_BAR = new KeyValuePair().withName('BAR').withValue('2')
        def task = Mock(TaskRun)
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig(memory: '8GB', cpus: 4, maxRetries: 2)

        def handler = Spy(AwsBatchTaskHandler)

        when:
        def req = handler.newSubmitRequest(task)
        then:
        1 * handler.getAwsOptions() >> { new AwsOptions(cliPath: '/bin/aws') }
        1 * handler.getJobQueue(task) >> 'queue1'
        1 * handler.getJobDefinition(task) >> 'job-def:1'
        1 * handler.getEnvironmentVars() >> [VAR_FOO, VAR_BAR]
        1 * handler.wrapperFile >> Paths.get('/bucket/test/.command.run')
        1 * handler.getLogFile() >> Paths.get('/bucket/test/.command.log')

        req.getJobName() == 'batchtask'
        req.getJobQueue() == 'queue1'
        req.getJobDefinition() == 'job-def:1'
        // no error `retry` error strategy is defined by NF, use `maxRetries` to se Batch attempts
        req.getRetryStrategy() == new RetryStrategy().withAttempts(3)
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
        1 * handler.resolveJobDefinition(IMAGE) >> JOB_NAME
        result == JOB_NAME

    }

    def 'should return job envs'() {
        given:
        def VAR_FOO = new KeyValuePair().withName('FOO').withValue('hello')
        def VAR_BAR = new KeyValuePair().withName('BAR').withValue('world')
        def VAR_NXF = new KeyValuePair().withName('NXF_DEBUG').withValue('2')

        def bean = new TaskBean()
        bean.environment = [FOO:'hello', BAR: 'world']
        def handler = [:] as AwsBatchTaskHandler
        handler.bean = bean

        when:
        def vars = handler.getEnvironmentVars()
        then:
        vars.size() == 0

        when:
        handler.environment = [NXF_DEBUG: 2, NXF_FOO: 'ignore']
        vars = handler.getEnvironmentVars()
        then:
        vars == [ VAR_NXF ]
    }

    def 'should strip invalid chars for job definition name' () {
        given:
        def handler = Spy(AwsBatchTaskHandler)

        expect:
        handler.normalizeJobDefinitionName(null) == null
        handler.normalizeJobDefinitionName('foo') == 'nf-foo'
        handler.normalizeJobDefinitionName('foo:1') == 'nf-foo-1'
    }


    def 'should validate job validation method' () {
        given:
        def IMAGE = 'foo/bar:1.0'
        def JOB_NAME = 'nf-foo-bar-1-0'
        def JOB_ID= '123'
        def handler = Spy(AwsBatchTaskHandler)

        def req = Mock(RegisterJobDefinitionRequest)

        when:
        handler.resolveJobDefinition(IMAGE)
        then:
        1 * handler.makeJobDefRequest(IMAGE) >> req
        1 * req.getJobDefinitionName() >> JOB_NAME
        2 * req.getParameters() >> [ 'nf-token': JOB_ID ]
        1 * handler.findJobDef(JOB_NAME, JOB_ID) >> null
        1 * handler.createJobDef(req) >> null

        when:
        handler.resolveJobDefinition(IMAGE)
        then:
        // second time are not invoked for the same image
        0 * handler.makeJobDefRequest(IMAGE) >> req
        0 * handler.findJobDef(JOB_NAME, JOB_ID) >> null
        0 * handler.createJobDef(req) >> null

    }

    def 'should verify job definition existence' () {

        given:
        def JOB_NAME = 'foo-bar-1-0'
        def JOB_ID = '123'
        def client = Mock(AWSBatchClient)
        def handler = Spy(AwsBatchTaskHandler)
        handler.client = client

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
        2 * job.getParameters() >> ['nf-token': JOB_ID]
        1 * job.getRevision() >> 3
        result == "$JOB_NAME:3"

        when:
        result = handler.findJobDef(JOB_NAME, JOB_ID)
        then:
        1 * client.describeJobDefinitions(req) >> res
        1 * res.getJobDefinitions() >> [job]
        1 * job.getStatus() >> 'ACTIVE'
        2 * job.getParameters() >> [:]
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
        def client = Mock(AWSBatchClient)
        def handler = Spy(AwsBatchTaskHandler)
        handler.client = client

        def req = Mock(RegisterJobDefinitionRequest)
        def res = Mock(RegisterJobDefinitionResult)

        when:
        def result = handler.createJobDef(req)
        then:
        1 * client.registerJobDefinition(req) >> res
        1 * res.getJobDefinitionName() >> JOB_NAME
        1 * res.getRevision() >> 10
        result == "$JOB_NAME:10"

    }


    def 'should create a job definition request object' () {
        given:
        def IMAGE = 'foo/bar:1.0'
        def JOB_NAME = 'nf-foo-bar-1-0'
        def handler = Spy(AwsBatchTaskHandler)

        when:
        def result = handler.makeJobDefRequest(IMAGE)
        then:
        1 * handler.normalizeJobDefinitionName(IMAGE) >> JOB_NAME
        1 * handler.getAwsOptions() >> new AwsOptions()
        result.jobDefinitionName == JOB_NAME
        result.type == 'container'
        result.parameters.'nf-token' == 'fdb5ef295f566138a43252b2ea272282'
        !result.containerProperties.mountPoints

        when:
        result = handler.makeJobDefRequest(IMAGE)
        then:
        1 * handler.normalizeJobDefinitionName(IMAGE) >> JOB_NAME
        1 * handler.getAwsOptions() >> new AwsOptions(cliPath: '/home/conda/bin/aws')
        result.jobDefinitionName == JOB_NAME
        result.type == 'container'
        result.parameters.'nf-token' == '9c56fd073d32e0c29f51f12afdfe4750'
        result.containerProperties.mountPoints[0].sourceVolume == 'aws-cli'
        result.containerProperties.mountPoints[0].containerPath == '/home/conda'
        result.containerProperties.mountPoints[0].readOnly
        result.containerProperties.volumes[0].host.sourcePath == '/home/conda'
        result.containerProperties.volumes[0].name == 'aws-cli'

    }


}
