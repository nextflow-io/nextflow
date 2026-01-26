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
import io.seqera.sched.api.schema.v1a1.PriceModel as SchedPriceModel
import io.seqera.sched.api.schema.v1a1.ProvisioningModel
import nextflow.cloud.types.PriceModel
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

}
