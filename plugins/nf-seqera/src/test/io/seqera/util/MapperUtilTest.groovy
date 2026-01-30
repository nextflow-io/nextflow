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
 *
 */

package io.seqera.util

import io.seqera.config.MachineRequirementOpts
import io.seqera.sched.api.schema.v1a1.DiskAllocation
import io.seqera.sched.api.schema.v1a1.PriceModel as SchedPriceModel
import io.seqera.sched.api.schema.v1a1.ProvisioningModel
import nextflow.cloud.types.PriceModel
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
        MapperUtil.toMachineRequirement(null) == null
    }

    def 'should return null for empty opts' () {
        expect:
        MapperUtil.toMachineRequirement(new MachineRequirementOpts([:])) == null
    }

    def 'should map arch only' () {
        when:
        def result = MapperUtil.toMachineRequirement(new MachineRequirementOpts([arch: 'arm64']))

        then:
        result.arch == 'arm64'
        result.provisioning == null
        result.maxSpotAttempts == null
        result.machineFamilies == null
    }

    def 'should map all fields' () {
        when:
        def result = MapperUtil.toMachineRequirement(new MachineRequirementOpts([
            arch: 'x86_64',
            provisioning: 'spotFirst',
            maxSpotAttempts: 3,
            machineFamilies: ['m5', 'c5']
        ]))

        then:
        result.arch == 'x86_64'
        result.provisioning == ProvisioningModel.SPOT_FIRST
        result.maxSpotAttempts == 3
        result.machineFamilies == ['m5', 'c5']
    }

    def 'should map provisioning model' () {
        expect:
        MapperUtil.toProvisioningModel(null) == null
        MapperUtil.toProvisioningModel('spot') == ProvisioningModel.SPOT
        MapperUtil.toProvisioningModel('ondemand') == ProvisioningModel.ONDEMAND
        MapperUtil.toProvisioningModel('spotFirst') == ProvisioningModel.SPOT_FIRST
    }

    def 'should map price model' () {
        expect:
        MapperUtil.toPriceModel(null) == null
        MapperUtil.toPriceModel(SchedPriceModel.SPOT) == PriceModel.spot
        MapperUtil.toPriceModel(SchedPriceModel.STANDARD) == PriceModel.standard
    }

    // tests for toMachineRequirement with task arch

    def 'should return null when both opts and taskArch are null' () {
        expect:
        MapperUtil.toMachineRequirement(null, null) == null
    }

    def 'should use taskArch when opts is null' () {
        when:
        def result = MapperUtil.toMachineRequirement(null, 'arm64')

        then:
        result.arch == 'arm64'
        result.provisioning == null
    }

    def 'should use taskArch over config arch' () {
        when:
        def result = MapperUtil.toMachineRequirement(
            new MachineRequirementOpts([arch: 'x86_64', provisioning: 'spot']),
            'arm64'
        )

        then:
        result.arch == 'arm64'
        result.provisioning == ProvisioningModel.SPOT
    }

    def 'should use config arch when taskArch is null' () {
        when:
        def result = MapperUtil.toMachineRequirement(
            new MachineRequirementOpts([arch: 'x86_64', provisioning: 'spot']),
            null
        )

        then:
        result.arch == 'x86_64'
        result.provisioning == ProvisioningModel.SPOT
    }

    def 'should merge config settings with taskArch' () {
        when:
        def result = MapperUtil.toMachineRequirement(
            new MachineRequirementOpts([
                provisioning: 'spotFirst',
                maxSpotAttempts: 3,
                machineFamilies: ['m5', 'c5']
            ]),
            'arm64'
        )

        then:
        result.arch == 'arm64'
        result.provisioning == ProvisioningModel.SPOT_FIRST
        result.maxSpotAttempts == 3
        result.machineFamilies == ['m5', 'c5']
    }

    // tests for disk requirement mapping

    def 'should return null disk requirement for null disk size' () {
        expect:
        MapperUtil.toDiskRequirement(null) == null
    }

    def 'should return null disk requirement for zero disk size' () {
        expect:
        MapperUtil.toDiskRequirement(MemoryUnit.of(0)) == null
    }

    def 'should map disk size to disk requirement with defaults' () {
        when:
        def result = MapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'))

        then: 'disk size is set'
        result.sizeGiB == 100
        and: 'allocation defaults to null (task allocation on API side)'
        result.allocation == null
        and: 'EBS defaults are applied for task allocation'
        result.volumeType == MapperUtil.DEFAULT_DISK_TYPE
        result.throughputMiBps == MapperUtil.DEFAULT_DISK_THROUGHPUT_MIBPS
        result.encrypted == false
        result.iops == null
    }

    def 'should map disk size in different units' () {
        expect:
        MapperUtil.toDiskRequirement(MemoryUnit.of('1 TB')).sizeGiB == 1024
        MapperUtil.toDiskRequirement(MemoryUnit.of('50 GB')).sizeGiB == 50
        and: 'defaults are applied'
        MapperUtil.toDiskRequirement(MemoryUnit.of('1 TB')).volumeType == MapperUtil.DEFAULT_DISK_TYPE
        MapperUtil.toDiskRequirement(MemoryUnit.of('1 TB')).throughputMiBps == MapperUtil.DEFAULT_DISK_THROUGHPUT_MIBPS
    }

    def 'should include disk in machine requirement' () {
        when:
        def result = MapperUtil.toMachineRequirement(
            new MachineRequirementOpts([arch: 'x86_64']),
            null,
            MemoryUnit.of('200 GB')
        )

        then:
        result.arch == 'x86_64'
        result.disk != null
        result.disk.sizeGiB == 200
    }

    def 'should return machine requirement with only disk' () {
        when:
        def result = MapperUtil.toMachineRequirement(null, null, MemoryUnit.of('100 GB'))

        then:
        result != null
        result.arch == null
        result.disk != null
        result.disk.sizeGiB == 100
    }

    def 'should return null when no arch, no opts, and no disk' () {
        expect:
        MapperUtil.toMachineRequirement(null, null, null) == null
    }

    // tests for custom disk configuration options

    def 'should throw exception for invalid disk type' () {
        given:
        def opts = new MachineRequirementOpts([diskType: 'local/nvme'])

        when:
        MapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains("Invalid disk type: local/nvme")
        e.message.contains("Supported types:")
    }

    def 'should use custom disk type from config' () {
        given:
        def opts = new MachineRequirementOpts([diskType: 'ebs/io1', diskIops: 10000])

        when:
        def result = MapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        result.sizeGiB == 100
        result.volumeType == 'ebs/io1'
        result.iops == 10000
        result.throughputMiBps == null  // throughput only for gp3
    }

    def 'should use custom throughput from config' () {
        given:
        def opts = new MachineRequirementOpts([diskThroughputMiBps: 500])

        when:
        def result = MapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        result.volumeType == MapperUtil.DEFAULT_DISK_TYPE
        result.throughputMiBps == 500
    }

    def 'should use encryption from config' () {
        given:
        def opts = new MachineRequirementOpts([diskEncrypted: true])

        when:
        def result = MapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        result.encrypted == true
    }

    def 'should use all disk options from config' () {
        given:
        def opts = new MachineRequirementOpts([
            diskType: 'ebs/gp3',
            diskThroughputMiBps: 600,
            diskIops: 8000,
            diskEncrypted: true
        ])

        when:
        def result = MapperUtil.toDiskRequirement(MemoryUnit.of('200 GB'), opts)

        then:
        result.sizeGiB == 200
        result.volumeType == 'ebs/gp3'
        result.throughputMiBps == 600
        result.iops == 8000
        result.encrypted == true
    }

    def 'should pass disk options through machine requirement' () {
        given:
        def opts = new MachineRequirementOpts([
            arch: 'arm64',
            diskType: 'ebs/io2',
            diskIops: 15000,
            diskEncrypted: true
        ])

        when:
        def result = MapperUtil.toMachineRequirement(opts, null, MemoryUnit.of('500 GB'))

        then:
        result.arch == 'arm64'
        result.disk.sizeGiB == 500
        result.disk.volumeType == 'ebs/io2'
        result.disk.iops == 15000
        result.disk.encrypted == true
        result.disk.throughputMiBps == null  // io2 doesn't use throughput
    }

    // tests for disk allocation mapping

    def 'should map disk allocation' () {
        expect:
        MapperUtil.toDiskAllocation(null) == null
        MapperUtil.toDiskAllocation('task') == DiskAllocation.TASK
        MapperUtil.toDiskAllocation('node') == DiskAllocation.NODE
    }

    def 'should throw exception for invalid disk allocation' () {
        when:
        MapperUtil.toDiskAllocation('invalid')

        then:
        thrown(IllegalArgumentException)
    }

    def 'should use task disk allocation from config' () {
        given:
        def opts = new MachineRequirementOpts([diskAllocation: 'task'])

        when:
        def result = MapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        result.allocation == DiskAllocation.TASK
    }

    def 'should use node disk allocation from config' () {
        given:
        def opts = new MachineRequirementOpts([diskAllocation: 'node'])

        when:
        def result = MapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        result.allocation == DiskAllocation.NODE
        result.sizeGiB == 100
        result.volumeType == null  // node allocation doesn't set EBS options
        result.throughputMiBps == null
        result.iops == null
        result.encrypted == null
    }

    def 'should default to null disk allocation when not specified' () {
        when:
        def result = MapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'))

        then:
        result.allocation == null
    }

    def 'should include disk allocation in machine requirement' () {
        given:
        def opts = new MachineRequirementOpts([
            arch: 'x86_64',
            diskAllocation: 'node'
        ])

        when:
        def result = MapperUtil.toMachineRequirement(opts, null, MemoryUnit.of('200 GB'))

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
        MapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

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
        MapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

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
        MapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

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
        MapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

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
        MapperUtil.toDiskRequirement(MemoryUnit.of('100 GB'), opts)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('diskType')
        e.message.contains('diskIops')
        e.message.contains('diskEncrypted')
    }

}
