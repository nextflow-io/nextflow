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

package nextflow.cloud.aws.batch.model

import software.amazon.awssdk.services.batch.model.ContainerProperties
import software.amazon.awssdk.services.batch.model.EphemeralStorage
import software.amazon.awssdk.services.batch.model.KeyValuePair
import software.amazon.awssdk.services.batch.model.LinuxParameters
import software.amazon.awssdk.services.batch.model.LogConfiguration
import software.amazon.awssdk.services.batch.model.MountPoint
import software.amazon.awssdk.services.batch.model.NetworkConfiguration
import software.amazon.awssdk.services.batch.model.ResourceRequirement
import software.amazon.awssdk.services.batch.model.ResourceType
import software.amazon.awssdk.services.batch.model.RuntimePlatform
import software.amazon.awssdk.services.batch.model.Ulimit
import software.amazon.awssdk.services.batch.model.Volume
import spock.lang.Specification

/**
 * @author Nextflow Authors
 */
class ContainerPropertiesModelTest extends Specification {

    def 'should create empty model'() {
        when:
        def model = new ContainerPropertiesModel()

        then:
        model.image == null
        model.command == null
        model.resourceRequirements == null
        model.jobRoleArn == null
        model.executionRoleArn == null
        model.linuxParameters == null
        model.environment == null
        model.privileged == false
        model.user == null
        model.readonlyRootFilesystem == false
        model.ulimits == null
        model.logConfiguration == null
        model.mountPoints == null
        model.volumes == null
        model.networkConfiguration == null
        model.ephemeralStorage == null
        model.runtimePlatform == null
    }

    def 'should set and get image'() {
        given:
        def model = new ContainerPropertiesModel()

        when:
        def result = model.image('ubuntu:20.04')

        then:
        result == model
        model.image == 'ubuntu:20.04'
    }

    def 'should set and get command'() {
        given:
        def model = new ContainerPropertiesModel()

        when:
        def result = model.command('echo', 'hello', 'world')

        then:
        result == model
        model.command == ['echo', 'hello', 'world']
        model.command.size() == 3
    }

    def 'should set and get resource requirements'() {
        given:
        def model = new ContainerPropertiesModel()
        def req1 = ResourceRequirement.builder()
            .type(ResourceType.VCPU)
            .value('1')
            .build()
        def req2 = ResourceRequirement.builder()
            .type(ResourceType.MEMORY)
            .value('1024')
            .build()

        when:
        def result = model.resourceRequirements(req1, req2)

        then:
        result == model
        model.resourceRequirements.size() == 2
        model.resourceRequirements[0] == req1
        model.resourceRequirements[1] == req2
    }

    def 'should set and get job role arn'() {
        given:
        def model = new ContainerPropertiesModel()
        def arn = 'arn:aws:iam::123456789012:role/BatchJobRole'

        when:
        def result = model.jobRoleArn(arn)

        then:
        result == model
        model.jobRoleArn == arn
    }

    def 'should set and get execution role arn'() {
        given:
        def model = new ContainerPropertiesModel()
        def arn = 'arn:aws:iam::123456789012:role/BatchExecutionRole'

        when:
        def result = model.executionRoleArn(arn)

        then:
        result == model
        model.executionRoleArn == arn
    }

    def 'should set and get user'() {
        given:
        def model = new ContainerPropertiesModel()

        when:
        def result = model.user('batch-user')

        then:
        result == model
        model.user == 'batch-user'
    }

    def 'should set and get readonly root filesystem'() {
        given:
        def model = new ContainerPropertiesModel()

        when:
        def result = model.readonlyRootFilesystem(true)

        then:
        result == model
        model.readonlyRootFilesystem == true
    }

    def 'should set and get environment'() {
        given:
        def model = new ContainerPropertiesModel()
        def env = [
            KeyValuePair.builder().name('VAR1').value('value1').build(),
            KeyValuePair.builder().name('VAR2').value('value2').build()
        ] as ArrayList<KeyValuePair>

        when:
        def result = model.environment(env)

        then:
        result == model
        model.environment == env
        model.environment.size() == 2
        model.environment[0].name() == 'VAR1'
        model.environment[0].value() == 'value1'
        model.environment[1].name() == 'VAR2'
        model.environment[1].value() == 'value2'
    }

    def 'should set and get linux parameters'() {
        given:
        def model = new ContainerPropertiesModel()
        def linuxParams = LinuxParameters.builder()
            .initProcessEnabled(true)
            .build()

        when:
        def result = model.linuxParameters(linuxParams)

        then:
        result == model
        model.linuxParameters == linuxParams
    }

    def 'should set and get privileged'() {
        given:
        def model = new ContainerPropertiesModel()

        when:
        def result = model.privileged(true)

        then:
        result == model
        model.privileged == true
    }

    def 'should set and get ulimits'() {
        given:
        def model = new ContainerPropertiesModel()
        def ulimits = [
            Ulimit.builder().name('nofile').softLimit(1024).hardLimit(2048).build(),
            Ulimit.builder().name('nproc').softLimit(16).hardLimit(32).build()
        ] as ArrayList<Ulimit>

        when:
        def result = model.ulimits(ulimits)

        then:
        result == model
        model.ulimits == ulimits
        model.ulimits.size() == 2
        model.ulimits[0].name() == 'nofile'
        model.ulimits[0].softLimit() == 1024
        model.ulimits[0].hardLimit() == 2048
    }

    def 'should set and get log configuration'() {
        given:
        def model = new ContainerPropertiesModel()
        def logConfig = LogConfiguration.builder()
            .logDriver('awslogs')
            .options(['awslogs-group': '/aws/batch/job'])
            .build()

        when:
        def result = model.logConfiguration(logConfig)

        then:
        result == model
        model.logConfiguration == logConfig
    }

    def 'should set and get mount points'() {
        given:
        def model = new ContainerPropertiesModel()
        def mountPoints = [
            MountPoint.builder()
                .sourceVolume('tmp')
                .containerPath('/tmp')
                .readOnly(false)
                .build()
        ]

        when:
        def result = model.mountPoints(mountPoints)

        then:
        result == model
        model.mountPoints == mountPoints
        model.mountPoints.size() == 1
        model.mountPoints[0].sourceVolume() == 'tmp'
        model.mountPoints[0].containerPath() == '/tmp'
        model.mountPoints[0].readOnly() == false
    }

    def 'should set and get volumes'() {
        given:
        def model = new ContainerPropertiesModel()
        def volumes = [
            Volume.builder()
                .name('tmp')
                .build()
        ]

        when:
        def result = model.volumes(volumes)

        then:
        result == model
        model.volumes == volumes
        model.volumes.size() == 1
        model.volumes[0].name() == 'tmp'
    }

    def 'should set and get network configuration'() {
        given:
        def model = new ContainerPropertiesModel()
        def networkConfig = NetworkConfiguration.builder()
            .assignPublicIp('ENABLED')
            .build()

        when:
        def result = model.networkConfiguration(networkConfig)

        then:
        result == model
        model.networkConfiguration == networkConfig
    }

    def 'should set and get ephemeral storage'() {
        given:
        def model = new ContainerPropertiesModel()
        def ephemeralStorage = EphemeralStorage.builder()
            .sizeInGiB(20)
            .build()

        when:
        def result = model.ephemeralStorage(ephemeralStorage)

        then:
        result == model
        model.ephemeralStorage == ephemeralStorage
    }

    def 'should set and get runtime platform'() {
        given:
        def model = new ContainerPropertiesModel()
        def runtimePlatform = RuntimePlatform.builder()
            .operatingSystemFamily('LINUX')
            .cpuArchitecture('X86_64')
            .build()

        when:
        def result = model.runtimePlatform(runtimePlatform)

        then:
        result == model
        model.runtimePlatform == runtimePlatform
    }

    def 'should support method chaining'() {
        given:
        def model = new ContainerPropertiesModel()
        def req = ResourceRequirement.builder()
            .type(ResourceType.VCPU)
            .value('1')
            .build()
        def env = [
            KeyValuePair.builder().name('VAR1').value('value1').build()
        ] as ArrayList<KeyValuePair>

        when:
        def result = model
            .image('ubuntu:20.04')
            .command('echo', 'hello')
            .resourceRequirements(req)
            .jobRoleArn('arn:aws:iam::123456789012:role/BatchJobRole')
            .executionRoleArn('arn:aws:iam::123456789012:role/BatchExecutionRole')
            .user('batch-user')
            .readonlyRootFilesystem(true)
            .environment(env)
            .privileged(false)

        then:
        result == model
        model.image == 'ubuntu:20.04'
        model.command == ['echo', 'hello']
        model.resourceRequirements.size() == 1
        model.jobRoleArn == 'arn:aws:iam::123456789012:role/BatchJobRole'
        model.executionRoleArn == 'arn:aws:iam::123456789012:role/BatchExecutionRole'
        model.user == 'batch-user'
        model.readonlyRootFilesystem == true
        model.environment.size() == 1
        model.privileged == false
    }

    def 'should generate proper toString'() {
        given:
        def model = new ContainerPropertiesModel()
        def req = ResourceRequirement.builder()
            .type(ResourceType.VCPU)
            .value('1')
            .build()

        when:
        model.image('ubuntu:20.04')
             .command('echo', 'hello')
             .resourceRequirements(req)
             .jobRoleArn('arn:aws:iam::123456789012:role/BatchJobRole')
             .privileged(true)
             .user('batch-user')

        then:
        def toString = model.toString()
        toString.contains('ContainerPropertiesModel{')
        toString.contains("image='ubuntu:20.04'")
        toString.contains('command=[echo, hello]')
        toString.contains('resourceRequirements=')
        toString.contains("jobRoleArn='arn:aws:iam::123456789012:role/BatchJobRole'")
        toString.contains('privileged=true')
        toString.contains("user='batch-user'")
    }

    def 'should handle null values in toString'() {
        given:
        def model = new ContainerPropertiesModel()

        when:
        def toString = model.toString()

        then:
        toString.contains('ContainerPropertiesModel{')
        toString.contains("image='null'")
        toString.contains('command=null')
        toString.contains('resourceRequirements=null')
        toString.contains("jobRoleArn='null'")
        toString.contains('privileged=false')
        toString.contains("user='null'")
        toString.contains('readonlyRootFilesystem=false')
    }

    def 'should handle empty collections'() {
        given:
        def model = new ContainerPropertiesModel()

        when:
        model.environment([] as ArrayList<KeyValuePair>)
             .ulimits([] as ArrayList<Ulimit>)
             .mountPoints([])
             .volumes([])

        then:
        model.environment == []
        model.ulimits == []
        model.mountPoints == []
        model.volumes == []
    }

    def 'should handle single command argument'() {
        given:
        def model = new ContainerPropertiesModel()

        when:
        model.command('single-command')

        then:
        model.command == ['single-command']
        model.command.size() == 1
    }

    def 'should handle single resource requirement'() {
        given:
        def model = new ContainerPropertiesModel()
        def req = ResourceRequirement.builder()
            .type(ResourceType.MEMORY)
            .value('512')
            .build()

        when:
        model.resourceRequirements(req)

        then:
        model.resourceRequirements.size() == 1
        model.resourceRequirements[0] == req
    }

    def 'should handle boolean values correctly'() {
        given:
        def model = new ContainerPropertiesModel()

        when:
        model.privileged(false)
             .readonlyRootFilesystem(false)

        then:
        model.privileged == false
        model.readonlyRootFilesystem == false

        when:
        model.privileged(true)
             .readonlyRootFilesystem(true)

        then:
        model.privileged == true
        model.readonlyRootFilesystem == true
    }

    def 'should convert to ContainerProperties with all fields'() {
        given:
        def model = new ContainerPropertiesModel()
        def req = ResourceRequirement.builder()
            .type(ResourceType.VCPU)
            .value('1')
            .build()
        def env = [
            KeyValuePair.builder().name('VAR1').value('value1').build()
        ] as ArrayList<KeyValuePair>
        def ulimits = [
            Ulimit.builder().name('nofile').softLimit(1024).hardLimit(2048).build()
        ] as ArrayList<Ulimit>
        def logConfig = LogConfiguration.builder()
            .logDriver('awslogs')
            .build()
        def mountPoints = [
            MountPoint.builder()
                .sourceVolume('tmp')
                .containerPath('/tmp')
                .build()
        ]
        def volumes = [
            Volume.builder()
                .name('tmp')
                .build()
        ]
        def networkConfig = NetworkConfiguration.builder()
            .assignPublicIp('ENABLED')
            .build()
        def ephemeralStorage = EphemeralStorage.builder()
            .sizeInGiB(20)
            .build()
        def runtimePlatform = RuntimePlatform.builder()
            .operatingSystemFamily('LINUX')
            .build()
        def linuxParams = LinuxParameters.builder()
            .initProcessEnabled(true)
            .build()

        when:
        model.image('ubuntu:20.04')
             .command('echo', 'hello')
             .resourceRequirements(req)
             .jobRoleArn('arn:aws:iam::123456789012:role/BatchJobRole')
             .executionRoleArn('arn:aws:iam::123456789012:role/BatchExecutionRole')
             .linuxParameters(linuxParams)
             .environment(env)
             .privileged(true)
             .user('batch-user')
             .readonlyRootFilesystem(true)
             .ulimits(ulimits)
             .logConfiguration(logConfig)
             .mountPoints(mountPoints)
             .volumes(volumes)
             .networkConfiguration(networkConfig)
             .ephemeralStorage(ephemeralStorage)
             .runtimePlatform(runtimePlatform)

        def containerProperties = model.toBatchContainerProperties()

        then:
        containerProperties instanceof ContainerProperties
        containerProperties.image() == 'ubuntu:20.04'
        containerProperties.command() == ['echo', 'hello']
        containerProperties.resourceRequirements().size() == 1
        containerProperties.resourceRequirements()[0] == req
        containerProperties.jobRoleArn() == 'arn:aws:iam::123456789012:role/BatchJobRole'
        containerProperties.executionRoleArn() == 'arn:aws:iam::123456789012:role/BatchExecutionRole'
        containerProperties.linuxParameters() == linuxParams
        containerProperties.environment().size() == 1
        containerProperties.environment()[0].name() == 'VAR1'
        containerProperties.privileged() == true
        containerProperties.user() == 'batch-user'
        containerProperties.readonlyRootFilesystem() == true
        containerProperties.ulimits().size() == 1
        containerProperties.ulimits()[0].name() == 'nofile'
        containerProperties.logConfiguration() == logConfig
        containerProperties.mountPoints().size() == 1
        containerProperties.mountPoints()[0].sourceVolume() == 'tmp'
        containerProperties.volumes().size() == 1
        containerProperties.volumes()[0].name() == 'tmp'
        containerProperties.networkConfiguration() == networkConfig
        containerProperties.ephemeralStorage() == ephemeralStorage
        containerProperties.runtimePlatform() == runtimePlatform
    }

    def 'should convert to ContainerProperties with null fields'() {
        given:
        def model = new ContainerPropertiesModel()

        when:
        def containerProperties = model.toBatchContainerProperties()

        then:
        containerProperties instanceof ContainerProperties
        containerProperties.image() == null
        containerProperties.command() == []
        containerProperties.resourceRequirements() == []
        containerProperties.jobRoleArn() == null
        containerProperties.executionRoleArn() == null
        containerProperties.linuxParameters() == null
        containerProperties.environment() == []
        containerProperties.privileged() == null
        containerProperties.user() == null
        containerProperties.readonlyRootFilesystem() == null
        containerProperties.ulimits() == []
        containerProperties.logConfiguration() == null
        containerProperties.mountPoints() == []
        containerProperties.volumes() == []
        containerProperties.networkConfiguration() == null
        containerProperties.ephemeralStorage() == null
        containerProperties.runtimePlatform() == null
    }

    def 'should convert to ContainerProperties with empty collections'() {
        given:
        def model = new ContainerPropertiesModel()

        when:
        model.environment([] as ArrayList<KeyValuePair>)
             .ulimits([] as ArrayList<Ulimit>)
             .mountPoints([])
             .volumes([])

        def containerProperties = model.toBatchContainerProperties()

        then:
        containerProperties instanceof ContainerProperties
        containerProperties.environment() == []
        containerProperties.ulimits() == []
        containerProperties.mountPoints() == []
        containerProperties.volumes() == []
    }
}
