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

package io.seqera.util

import io.seqera.config.MachineRequirementOpts
import io.seqera.sched.api.schema.v1a1.DiskAllocation
import io.seqera.sched.api.schema.v1a1.EcsCapacityMode
import io.seqera.sched.api.schema.v1a1.PriceModel as SchedPriceModel
import io.seqera.sched.api.schema.v1a1.ProvisioningModel
import nextflow.cloud.types.PriceModel
import nextflow.fusion.FusionConfig
import nextflow.util.MemoryUnit
import spock.lang.Specification

/**
 * Unit tests for SchemaMapper
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class MapperUtilTest extends Specification {

    def 'should return null for null opts' () {
        expect:
        SchemaMapperUtil.toMachineRequirement(null) == null
    }

    def 'should return null for empty opts' () {
        expect:
        SchemaMapperUtil.toMachineRequirement(new MachineRequirementOpts([:])) == null
    }

    def 'should map arch only' () {
        when:
        def result = SchemaMapperUtil.toMachineRequirement(new MachineRequirementOpts([arch: 'arm64']))

        then:
        result.arch == 'arm64'
        result.provisioning == null
        result.maxSpotAttempts == null
        result.machineTypes == null
    }

    def 'should map all fields' () {
        when:
        def result = SchemaMapperUtil.toMachineRequirement(new MachineRequirementOpts([
            arch: 'x86_64',
            provisioning: 'spotFirst',
            maxSpotAttempts: 3,
            machineTypes: ['m5', 'c5']
        ]))

        then:
        result.arch == 'x86_64'
        result.provisioning == ProvisioningModel.SPOT_FIRST
        result.maxSpotAttempts == 3
        result.machineTypes == ['m5', 'c5']
    }

    def 'should map provisioning model' () {
        expect:
        SchemaMapperUtil.toProvisioningModel(null) == null
        SchemaMapperUtil.toProvisioningModel('spot') == ProvisioningModel.SPOT
        SchemaMapperUtil.toProvisioningModel('ondemand') == ProvisioningModel.ONDEMAND
        SchemaMapperUtil.toProvisioningModel('spotFirst') == ProvisioningModel.SPOT_FIRST
    }

    def 'should map price model' () {
        expect:
        SchemaMapperUtil.toPriceModel(null) == null
        SchemaMapperUtil.toPriceModel(SchedPriceModel.SPOT) == PriceModel.spot
        SchemaMapperUtil.toPriceModel(SchedPriceModel.STANDARD) == PriceModel.standard
    }

    // tests for toMachineRequirement with task arch

    def 'should return null when both opts and taskArch are null' () {
        expect:
        SchemaMapperUtil.toMachineRequirement(null, null) == null
    }

    def 'should use taskArch when opts is null' () {
        when:
        def result = SchemaMapperUtil.toMachineRequirement(null, 'arm64')

        then:
        result.arch == 'arm64'
        result.provisioning == null
    }

    def 'should use taskArch over config arch' () {
        when:
        def result = SchemaMapperUtil.toMachineRequirement(
            new MachineRequirementOpts([arch: 'x86_64', provisioning: 'spot']),
            'arm64'
        )

        then:
        result.arch == 'arm64'
        result.provisioning == ProvisioningModel.SPOT
    }

    def 'should use config arch when taskArch is null' () {
        when:
        def result = SchemaMapperUtil.toMachineRequirement(
            new MachineRequirementOpts([arch: 'x86_64', provisioning: 'spot']),
            null
        )

        then:
        result.arch == 'x86_64'
        result.provisioning == ProvisioningModel.SPOT
    }

    def 'should merge config settings with taskArch' () {
        when:
        def result = SchemaMapperUtil.toMachineRequirement(
            new MachineRequirementOpts([
                provisioning: 'spotFirst',
                maxSpotAttempts: 3,
                machineTypes: ['m5', 'c5']
            ]),
            'arm64'
        )

        then:
        result.arch == 'arm64'
        result.provisioning == ProvisioningModel.SPOT_FIRST
        result.maxSpotAttempts == 3
        result.machineTypes == ['m5', 'c5']
    }

    // tests for disk requirement mapping

    def 'should return null disk requirement for null disk size' () {
        expect:
        SchemaMapperUtil.toDiskRequirement(null) == null
    }

    def 'should return null disk requirement for zero disk size' () {
        expect:
        SchemaMapperUtil.toDiskRequirement(MemoryUnit.of(0)) == null
    }

    def 'should map disk size to disk requirement with defaults' () {
        when:
        def result = SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'))

        then: 'disk size is set'
        result.sizeGiB == 100
        and: 'allocation defaults to node'
        result.allocation == DiskAllocation.NODE
        and: 'node allocation does not set EBS options'
        result.volumeType == null
        result.throughputMiBps == null
        result.encrypted == null
        result.iops == null
    }

    def 'should map disk size in different units' () {
        expect:
        SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('1 TB')).sizeGiB == 1024
        SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('50 GB')).sizeGiB == 50
        and: 'defaults to node allocation'
        SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('1 TB')).allocation == DiskAllocation.NODE
    }

    def 'should include disk in machine requirement' () {
        when:
        def result = SchemaMapperUtil.toMachineRequirement(
            new MachineRequirementOpts([arch: 'x86_64']),
            null,
            MemoryUnit.of('200 GB'),
            false
        )

        then:
        result.arch == 'x86_64'
        result.disk != null
        result.disk.sizeGiB == 200
    }

    def 'should return machine requirement with only disk' () {
        when:
        def result = SchemaMapperUtil.toMachineRequirement(null, null, MemoryUnit.of('100 GB'), false)

        then:
        result != null
        result.arch == null
        result.disk != null
        result.disk.sizeGiB == 100
    }

    def 'should return null when no arch, no opts, and no disk' () {
        expect:
        SchemaMapperUtil.toMachineRequirement(null, null, null, false) == null
    }

    // tests for custom disk configuration options

    def 'should throw exception for invalid disk type' () {
        given:
        def opts = new MachineRequirementOpts([diskAllocation: 'task', diskType: 'local/nvme'])

        when:
        SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains("Invalid disk type: local/nvme")
        e.message.contains("Supported types:")
    }

    def 'should use custom disk type from config' () {
        given:
        def opts = new MachineRequirementOpts([diskAllocation: 'task', diskType: 'ebs/io1', diskIops: 10000])

        when:
        def result = SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        result.sizeGiB == 100
        result.allocation == DiskAllocation.TASK
        result.volumeType == 'ebs/io1'
        result.iops == 10000
        result.throughputMiBps == null  // throughput only for gp3
    }

    def 'should use custom throughput from config' () {
        given:
        def opts = new MachineRequirementOpts([diskAllocation: 'task', diskThroughputMiBps: 500])

        when:
        def result = SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        result.allocation == DiskAllocation.TASK
        result.volumeType == SchemaMapperUtil.DEFAULT_DISK_TYPE
        result.throughputMiBps == 500
    }

    def 'should use encryption from config' () {
        given:
        def opts = new MachineRequirementOpts([diskAllocation: 'task', diskEncrypted: true])

        when:
        def result = SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        result.allocation == DiskAllocation.TASK
        result.encrypted == true
    }

    def 'should use all disk options from config' () {
        given:
        def opts = new MachineRequirementOpts([
            diskAllocation: 'task',
            diskType: 'ebs/gp3',
            diskThroughputMiBps: 600,
            diskIops: 8000,
            diskEncrypted: true
        ])

        when:
        def result = SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('200 GB'), opts)

        then:
        result.sizeGiB == 200
        result.allocation == DiskAllocation.TASK
        result.volumeType == 'ebs/gp3'
        result.throughputMiBps == 600
        result.iops == 8000
        result.encrypted == true
    }

    def 'should pass disk options through machine requirement' () {
        given:
        def opts = new MachineRequirementOpts([
            arch: 'arm64',
            diskAllocation: 'task',
            diskType: 'ebs/io2',
            diskIops: 15000,
            diskEncrypted: true
        ])

        when:
        def result = SchemaMapperUtil.toMachineRequirement(opts, null, MemoryUnit.of('500 GB'), false)

        then:
        result.arch == 'arm64'
        result.disk.sizeGiB == 500
        result.disk.allocation == DiskAllocation.TASK
        result.disk.volumeType == 'ebs/io2'
        result.disk.iops == 15000
        result.disk.encrypted == true
        result.disk.throughputMiBps == null  // io2 doesn't use throughput
    }

    // tests for disk allocation mapping

    def 'should map disk allocation' () {
        expect:
        SchemaMapperUtil.toDiskAllocation(null) == null
        SchemaMapperUtil.toDiskAllocation('task') == DiskAllocation.TASK
        SchemaMapperUtil.toDiskAllocation('node') == DiskAllocation.NODE
    }

    def 'should throw exception for invalid disk allocation' () {
        when:
        SchemaMapperUtil.toDiskAllocation('invalid')

        then:
        thrown(IllegalArgumentException)
    }

    def 'should use task disk allocation from config' () {
        given:
        def opts = new MachineRequirementOpts([diskAllocation: 'task'])

        when:
        def result = SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        result.allocation == DiskAllocation.TASK
    }

    def 'should use node disk allocation from config' () {
        given:
        def opts = new MachineRequirementOpts([diskAllocation: 'node'])

        when:
        def result = SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        result.allocation == DiskAllocation.NODE
        result.sizeGiB == 100
        result.volumeType == null  // node allocation doesn't set EBS options
        result.throughputMiBps == null
        result.iops == null
        result.encrypted == null
    }

    def 'should default to node disk allocation when not specified' () {
        when:
        def result = SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'))

        then:
        result.allocation == DiskAllocation.NODE
    }

    def 'should include disk allocation in machine requirement' () {
        given:
        def opts = new MachineRequirementOpts([
            arch: 'x86_64',
            diskAllocation: 'node'
        ])

        when:
        def result = SchemaMapperUtil.toMachineRequirement(opts, null, MemoryUnit.of('200 GB'), false)

        then:
        result.arch == 'x86_64'
        result.disk.sizeGiB == 200
        result.disk.allocation == DiskAllocation.NODE
    }

    // tests for node allocation validation

    def 'should throw error when diskType is set with node allocation' () {
        given:
        def opts = new MachineRequirementOpts([
            diskAllocation: 'node',
            diskType: 'ebs/gp3'
        ])

        when:
        SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('diskType')
        e.message.contains('node')
    }

    def 'should throw error when diskIops is set with node allocation' () {
        given:
        def opts = new MachineRequirementOpts([
            diskAllocation: 'node',
            diskIops: 10000
        ])

        when:
        SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('diskIops')
    }

    def 'should throw error when diskThroughputMiBps is set with node allocation' () {
        given:
        def opts = new MachineRequirementOpts([
            diskAllocation: 'node',
            diskThroughputMiBps: 500
        ])

        when:
        SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('diskThroughputMiBps')
    }

    def 'should throw error when diskEncrypted is set with node allocation' () {
        given:
        def opts = new MachineRequirementOpts([
            diskAllocation: 'node',
            diskEncrypted: true
        ])

        when:
        SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('diskEncrypted')
    }

    def 'should report all invalid options with node allocation' () {
        given:
        def opts = new MachineRequirementOpts([
            diskAllocation: 'node',
            diskType: 'ebs/io1',
            diskIops: 10000,
            diskEncrypted: true
        ])

        when:
        SchemaMapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('diskType')
        e.message.contains('diskIops')
        e.message.contains('diskEncrypted')
    }

    // tests for snapshot maxSpotAttempts defaulting

    def 'should return machine requirement with only snapshot enabled' () {
        when:
        def result = SchemaMapperUtil.toMachineRequirement(null, null, null, true)

        then:
        result != null
        result.snapshotEnabled == true
        result.maxSpotAttempts == FusionConfig.DEFAULT_SNAPSHOT_MAX_SPOT_ATTEMPTS
    }

    def 'should use explicit maxSpotAttempts when snapshot enabled' () {
        when:
        def result = SchemaMapperUtil.toMachineRequirement(new MachineRequirementOpts([maxSpotAttempts: 2]), null, null, true)

        then:
        result.snapshotEnabled == true
        result.maxSpotAttempts == 2
    }

    def 'should not default maxSpotAttempts when snapshot disabled' () {
        when:
        def result = SchemaMapperUtil.toMachineRequirement(new MachineRequirementOpts([arch: 'x86_64']), null, null, false)

        then:
        result.snapshotEnabled == null
        result.maxSpotAttempts == null
    }

    // tests for capacity mode mapping

    def 'should map capacity mode' () {
        expect:
        SchemaMapperUtil.toEcsCapacityMode(null) == null
        SchemaMapperUtil.toEcsCapacityMode('managed') == EcsCapacityMode.MANAGED
        SchemaMapperUtil.toEcsCapacityMode('asg') == EcsCapacityMode.ASG
    }

    def 'should throw exception for invalid capacity mode' () {
        when:
        SchemaMapperUtil.toEcsCapacityMode('invalid')

        then:
        thrown(IllegalArgumentException)
    }

    def 'should include capacity mode in machine requirement' () {
        when:
        def result = SchemaMapperUtil.toMachineRequirement(new MachineRequirementOpts([capacityMode: 'asg']))

        then:
        result != null
        result.capacityMode == EcsCapacityMode.ASG
    }

    def 'should include capacity mode in machine requirement with task arch' () {
        when:
        def result = SchemaMapperUtil.toMachineRequirement(
            new MachineRequirementOpts([capacityMode: 'managed', arch: 'arm64']),
            null,
            null,
            false
        )

        then:
        result.arch == 'arm64'
        result.capacityMode == EcsCapacityMode.MANAGED
    }

    def 'should combine snapshot with other machine requirement settings' () {
        when:
        def result = SchemaMapperUtil.toMachineRequirement(
            new MachineRequirementOpts([arch: 'arm64', provisioning: 'spot']),
            null,
            MemoryUnit.of('100 GB'),
            true
        )

        then:
        result.arch == 'arm64'
        result.provisioning == ProvisioningModel.SPOT
        result.disk.sizeGiB == 100
        result.snapshotEnabled == true
        result.maxSpotAttempts == FusionConfig.DEFAULT_SNAPSHOT_MAX_SPOT_ATTEMPTS
    }

}
