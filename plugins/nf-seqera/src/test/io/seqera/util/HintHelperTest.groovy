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
import nextflow.util.MemoryUnit
import spock.lang.Specification

/**
 * Tests for {@link HintHelper}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class HintHelperTest extends Specification {

    def 'should return base opts when no hints'() {
        given:
        def base = new MachineRequirementOpts([arch: 'x86_64'])

        when:
        def result = HintHelper.overlayHints(base, [:])
        then:
        result.arch == 'x86_64'
    }

    def 'should return base opts when no seqera hints'() {
        given:
        def base = new MachineRequirementOpts([arch: 'x86_64'])

        when:
        def result = HintHelper.overlayHints(base, [consumableResources: 'my-license'])
        then:
        result.arch == 'x86_64'
    }

    def 'should overlay arch hint'() {
        given:
        def base = new MachineRequirementOpts([arch: 'x86_64'])

        when:
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.arch': 'arm64'])
        then:
        result.arch == 'arm64'
    }

    def 'should overlay provisioning hint'() {
        given:
        def base = new MachineRequirementOpts([provisioning: 'ondemand'])

        when:
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.provisioning': 'spotFirst'])
        then:
        result.provisioning == 'spotFirst'
    }

    def 'should overlay maxSpotAttempts hint'() {
        given:
        def base = new MachineRequirementOpts([:])

        when:
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.maxSpotAttempts': 3])
        then:
        result.maxSpotAttempts == 3
    }

    def 'should overlay machineTypes from comma-separated string'() {
        given:
        def base = new MachineRequirementOpts([:])

        when:
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.machineTypes': 'm5,m5a,m6i'])
        then:
        result.machineTypes == ['m5', 'm5a', 'm6i']
    }

    def 'should overlay diskType hint'() {
        given:
        def base = new MachineRequirementOpts([:])

        when:
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.diskType': 'ebs/gp3'])
        then:
        result.diskType == 'ebs/gp3'
    }

    def 'should overlay diskThroughputMiBps hint'() {
        given:
        def base = new MachineRequirementOpts([:])

        when:
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.diskThroughputMiBps': 500])
        then:
        result.diskThroughputMiBps == 500
    }

    def 'should overlay diskIops hint'() {
        given:
        def base = new MachineRequirementOpts([:])

        when:
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.diskIops': 10000])
        then:
        result.diskIops == 10000
    }

    def 'should overlay diskEncrypted hint from string'() {
        given:
        def base = new MachineRequirementOpts([:])

        when:
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.diskEncrypted': 'true'])
        then:
        result.diskEncrypted == true
    }

    def 'should overlay diskAllocation hint'() {
        given:
        def base = new MachineRequirementOpts([:])

        when:
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.diskAllocation': 'node'])
        then:
        result.diskAllocation == 'node'
    }

    def 'should overlay diskMountPath hint'() {
        given:
        def base = new MachineRequirementOpts([:])

        when:
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.diskMountPath': '/data'])
        then:
        result.diskMountPath == '/data'
    }

    def 'should overlay diskSize hint'() {
        given:
        def base = new MachineRequirementOpts([:])

        when:
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.diskSize': '100.GB'])
        then:
        result.diskSize == MemoryUnit.of('100.GB')
    }

    def 'should overlay capacityMode hint'() {
        given:
        def base = new MachineRequirementOpts([:])

        when:
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.capacityMode': 'asg'])
        then:
        result.capacityMode == 'asg'
    }

    def 'should overlay multiple hints at once'() {
        given:
        def base = new MachineRequirementOpts([arch: 'x86_64', provisioning: 'ondemand'])

        when:
        def result = HintHelper.overlayHints(base, [
            'seqera/machineRequirement.arch': 'arm64',
            'seqera/machineRequirement.provisioning': 'spotFirst',
            'seqera/machineRequirement.maxSpotAttempts': 3,
            'seqera/machineRequirement.diskType': 'ebs/gp3'
        ])
        then:
        result.arch == 'arm64'
        result.provisioning == 'spotFirst'
        result.maxSpotAttempts == 3
        result.diskType == 'ebs/gp3'
    }

    def 'should preserve base values not overridden by hints'() {
        given:
        def base = new MachineRequirementOpts([arch: 'x86_64', provisioning: 'spot', diskType: 'ebs/gp3'])

        when:
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.arch': 'arm64'])
        then:
        result.arch == 'arm64'
        result.provisioning == 'spot'
        result.diskType == 'ebs/gp3'
    }

    def 'should error on unknown seqera hint'() {
        when:
        HintHelper.extractSeqeraHints(['seqera/machineRequirement.unknownField': 'value'])
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('Unknown Seqera Platform hint')
        e.message.contains('seqera/machineRequirement.unknownField')
    }

    def 'should ignore non-seqera hints in extraction'() {
        when:
        def result = HintHelper.extractSeqeraHints([
            consumableResources: 'my-license',
            'k8s/scheduling.nodeSelector': 'gpu=true',
            'seqera/machineRequirement.arch': 'arm64'
        ])
        then:
        result.size() == 1
        result['machineRequirement.arch'] == 'arm64'
    }

    def 'should handle null hints map'() {
        when:
        def result = HintHelper.extractSeqeraHints(null)
        then:
        result.isEmpty()
    }

}
