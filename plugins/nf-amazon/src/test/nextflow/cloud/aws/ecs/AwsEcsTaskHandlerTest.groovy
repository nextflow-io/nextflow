/*
 * Copyright 2020-2024, Seqera Labs
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
 *
 */

package nextflow.cloud.aws.ecs

import java.nio.file.Paths

import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.util.MemoryUnit
import software.amazon.awssdk.services.ecs.EcsClient
import software.amazon.awssdk.services.ecs.model.AssignPublicIp
import software.amazon.awssdk.services.ecs.model.Container
import software.amazon.awssdk.services.ecs.model.DescribeTasksResponse
import software.amazon.awssdk.services.ecs.model.Task
import spock.lang.Specification
import spock.lang.Unroll

/**
 * Tests for {@link AwsEcsTaskHandler}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsEcsTaskHandlerTest extends Specification {

    def 'should sanitize family name'() {
        given:
        def handler = [:] as AwsEcsTaskHandler

        expect:
        handler.sanitizeFamilyName('ubuntu:latest') == 'ubuntu-latest'
        handler.sanitizeFamilyName('quay.io/org/image:1.0') == 'quayio-org-image-10'
        handler.sanitizeFamilyName('registry.example.com/path/to/image@sha256:abc123') == 'registryexamplecom-path-to-image-sha256-abc123'
        handler.sanitizeFamilyName('simple') == 'simple'
        handler.sanitizeFamilyName(null) == 'unknown'
        handler.sanitizeFamilyName('') == 'unknown'
    }

    def 'should limit family name length'() {
        given:
        def handler = [:] as AwsEcsTaskHandler
        def longName = 'a' * 250

        when:
        def result = handler.sanitizeFamilyName(longName)

        then:
        result.size() == 200
    }

    @Unroll
    def 'should detect spot interruption: #stoppedReason'() {
        given:
        def handler = [:] as AwsEcsTaskHandler

        expect:
        handler.isSpotInterruption(stoppedReason, stopCode) == expected

        where:
        stoppedReason                                   | stopCode              | expected
        null                                            | null                  | false
        'Task finished successfully'                    | null                  | false
        'spot capacity unavailable'                     | null                  | true
        'Spot instance interruption'                    | null                  | true
        'Task stopped due to SPOT termination'          | null                  | true
        'Insufficient capacity in region'               | null                  | true
        'Normal termination'                            | 'SpotInterruption'    | true
        'User requested termination'                    | null                  | false
    }

    def 'should create task definition model'() {
        given:
        def task = Mock(TaskRun) {
            getContainer() >> 'ubuntu:latest'
            getConfig() >> new TaskConfig(cpus: 2, memory: new MemoryUnit('4 GB'))
            getWorkDir() >> Paths.get('/work/dir')
        }
        def awsOptions = Mock(AwsEcsOptions) {
            getExecutionRole() >> 'arn:aws:iam::123:role/exec'
            getTaskRole() >> 'arn:aws:iam::123:role/task'
            getLogsGroup() >> '/aws/ecs/test'
            getRegion() >> 'us-east-1'
        }
        def executor = Mock(AwsEcsExecutor) {
            getAwsOptions() >> awsOptions
        }
        def handler = new AwsEcsTaskHandler(task, executor)

        when:
        def model = handler.createTaskDefinitionModel()

        then:
        model.family == 'nf-ubuntu-latest'
        model.cpu == '2048'  // 2 CPUs * 1024
        model.memory == '4096'  // 4 GB in MiB
        model.executionRoleArn == 'arn:aws:iam::123:role/exec'
        model.taskRoleArn == 'arn:aws:iam::123:role/task'
        model.containerDefinitions.first().image == 'ubuntu:latest'
    }

    def 'should create task definition model with disk'() {
        given:
        def task = Mock(TaskRun) {
            getContainer() >> 'ubuntu:latest'
            getConfig() >> new TaskConfig(cpus: 1, memory: new MemoryUnit('2 GB'), disk: new MemoryUnit('100 GB'))
            getWorkDir() >> Paths.get('/work/dir')
        }
        def awsOptions = Mock(AwsEcsOptions) {
            getExecutionRole() >> 'arn:aws:iam::123:role/exec'
            getTaskRole() >> null
            getLogsGroup() >> '/aws/ecs/test'
            getRegion() >> 'us-east-1'
        }
        def executor = Mock(AwsEcsExecutor) {
            getAwsOptions() >> awsOptions
        }
        def handler = new AwsEcsTaskHandler(task, executor)

        when:
        def model = handler.createTaskDefinitionModel()

        then:
        model.ephemeralStorageGiB == 100
    }

    def 'should build run task request'() {
        given:
        def task = Mock(TaskRun) {
            getContainer() >> 'ubuntu:latest'
            getConfig() >> new TaskConfig(cpus: 2, memory: new MemoryUnit('4 GB'))
            getWorkDir() >> Paths.get('/work/dir')
        }
        def awsOptions = Mock(AwsEcsOptions) {
            getCluster() >> 'test-cluster'
            getAssignPublicIp() >> true
        }
        def executor = Mock(AwsEcsExecutor) {
            getAwsOptions() >> awsOptions
            getResolvedSubnets() >> ['subnet-123', 'subnet-456']
            getResolvedSecurityGroups() >> ['sg-789']
        }
        def handler = Spy(new AwsEcsTaskHandler(task, executor)) {
            fusionEnabled() >> false
            getEnvironment() >> [:]
            getContainerCommand() >> ['bash', '-c', 'echo hello']
        }
        handler.@taskDefArn = 'arn:aws:ecs:us-east-1:123:task-definition/test:1'

        when:
        def request = handler.newRunTaskRequest()

        then:
        request.cluster() == 'test-cluster'
        request.taskDefinition() == 'arn:aws:ecs:us-east-1:123:task-definition/test:1'
        request.count() == 1
        request.networkConfiguration().awsvpcConfiguration().subnets() == ['subnet-123', 'subnet-456']
        request.networkConfiguration().awsvpcConfiguration().securityGroups() == ['sg-789']
        request.networkConfiguration().awsvpcConfiguration().assignPublicIp() == AssignPublicIp.ENABLED
        request.overrides().containerOverrides().first().name() == 'main'
        request.overrides().containerOverrides().first().command() == ['bash', '-c', 'echo hello']
    }

    def 'should build run task request with disabled public IP'() {
        given:
        def task = Mock(TaskRun) {
            getContainer() >> 'ubuntu:latest'
            getConfig() >> new TaskConfig(cpus: 1, memory: new MemoryUnit('2 GB'))
            getWorkDir() >> Paths.get('/work/dir')
        }
        def awsOptions = Mock(AwsEcsOptions) {
            getCluster() >> 'test-cluster'
            getAssignPublicIp() >> false
        }
        def executor = Mock(AwsEcsExecutor) {
            getAwsOptions() >> awsOptions
            getResolvedSubnets() >> ['subnet-123']
            getResolvedSecurityGroups() >> ['sg-789']
        }
        def handler = Spy(new AwsEcsTaskHandler(task, executor)) {
            fusionEnabled() >> false
            getEnvironment() >> [:]
            getContainerCommand() >> ['bash', '-c', 'echo hello']
        }
        handler.@taskDefArn = 'arn:aws:ecs:us-east-1:123:task-definition/test:1'

        when:
        def request = handler.newRunTaskRequest()

        then:
        request.networkConfiguration().awsvpcConfiguration().assignPublicIp() == AssignPublicIp.DISABLED
    }

    def 'should check if running when task is running'() {
        given:
        def task = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/dir')
        }
        def ecsClient = Mock(EcsClient)
        def awsOptions = Mock(AwsEcsOptions) {
            getCluster() >> 'test-cluster'
        }
        def executor = Mock(AwsEcsExecutor) {
            getEcsClient() >> ecsClient
            getAwsOptions() >> awsOptions
        }
        def handler = new AwsEcsTaskHandler(task, executor)
        handler.@taskArn = 'arn:aws:ecs:us-east-1:123:task/test-cluster/abc123'
        handler.status = TaskStatus.SUBMITTED

        when:
        def result = handler.checkIfRunning()

        then:
        1 * ecsClient.describeTasks(_) >> DescribeTasksResponse.builder()
            .tasks(Task.builder()
                .taskArn('arn:aws:ecs:us-east-1:123:task/test-cluster/abc123')
                .lastStatus('RUNNING')
                .build())
            .build()

        and:
        result == true
        handler.status == TaskStatus.RUNNING
    }

    def 'should check if completed when task is stopped'() {
        given:
        def task = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/dir')
        }
        def ecsClient = Mock(EcsClient)
        def awsOptions = Mock(AwsEcsOptions) {
            getCluster() >> 'test-cluster'
            getMaxSpotAttempts() >> 5
        }
        def executor = Mock(AwsEcsExecutor) {
            getEcsClient() >> ecsClient
            getAwsOptions() >> awsOptions
        }
        def handler = new AwsEcsTaskHandler(task, executor)
        handler.@taskArn = 'arn:aws:ecs:us-east-1:123:task/test-cluster/abc123'
        handler.status = TaskStatus.RUNNING

        when:
        def result = handler.checkIfCompleted()

        then:
        1 * ecsClient.describeTasks(_) >> DescribeTasksResponse.builder()
            .tasks(Task.builder()
                .taskArn('arn:aws:ecs:us-east-1:123:task/test-cluster/abc123')
                .lastStatus('STOPPED')
                .containers(Container.builder()
                    .name('main')
                    .exitCode(0)
                    .build())
                .build())
            .build()

        and:
        result == true
        handler.status == TaskStatus.COMPLETED
    }

    def 'should return false when task not yet submitted'() {
        given:
        def task = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/dir')
        }
        def executor = Mock(AwsEcsExecutor)
        def handler = new AwsEcsTaskHandler(task, executor)
        // taskArn is null - not yet submitted

        when:
        def result = handler.checkIfRunning()

        then:
        result == false
    }

    def 'should return false when checkIfCompleted and not running'() {
        given:
        def task = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/dir')
        }
        def executor = Mock(AwsEcsExecutor)
        def handler = new AwsEcsTaskHandler(task, executor)
        handler.status = TaskStatus.SUBMITTED  // not yet RUNNING

        when:
        def result = handler.checkIfCompleted()

        then:
        result == false
    }

    def 'should enable Fusion FUSE support in task definition'() {
        given:
        def task = Mock(TaskRun) {
            getContainer() >> 'ubuntu:latest'
            getConfig() >> new TaskConfig(cpus: 1, memory: new MemoryUnit('2 GB'))
            getWorkDir() >> Paths.get('/work/dir')
        }
        def awsOptions = Mock(AwsEcsOptions) {
            getExecutionRole() >> 'arn:aws:iam::123:role/exec'
            getTaskRole() >> null
            getLogsGroup() >> '/aws/ecs/test'
            getRegion() >> 'us-east-1'
        }
        def executor = Mock(AwsEcsExecutor) {
            getAwsOptions() >> awsOptions
        }
        def handler = new AwsEcsTaskHandler(task, executor)

        when:
        def model = handler.createTaskDefinitionModel()

        then:
        // Verify fusionEnabled is set on container definition
        model.containerDefinitions.first().fusionEnabled == true
    }

    def 'should include Linux parameters for Fusion FUSE driver in container definition'() {
        given:
        def container = new nextflow.cloud.aws.ecs.model.ContainerDefinitionModel(name: 'main', image: 'ubuntu:latest')
        container.fusionEnabled = true

        when:
        def containerDef = container.toContainerDefinition()

        then:
        // Verify Linux parameters are set
        containerDef.linuxParameters() != null
        // Verify SYS_ADMIN capability is added
        containerDef.linuxParameters().capabilities().add().contains('SYS_ADMIN')
        // Verify /dev/fuse device is mounted
        containerDef.linuxParameters().devices().size() == 1
        containerDef.linuxParameters().devices().first().hostPath() == '/dev/fuse'
        containerDef.linuxParameters().devices().first().containerPath() == '/dev/fuse'
    }

    def 'should not include Linux parameters when Fusion not enabled'() {
        given:
        def container = new nextflow.cloud.aws.ecs.model.ContainerDefinitionModel(name: 'main', image: 'ubuntu:latest')
        container.fusionEnabled = false

        when:
        def containerDef = container.toContainerDefinition()

        then:
        // Verify no Linux parameters are set
        containerDef.linuxParameters() == null
    }
}
