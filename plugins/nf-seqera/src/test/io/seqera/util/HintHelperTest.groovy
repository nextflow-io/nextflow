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
        def base = new MachineRequirementOpts([provisioning: 'spot'])

        when:
        def result = HintHelper.overlayHints(base, [:])
        then:
        result.provisioning == 'spot'
    }

    def 'should return base opts when no seqera hints'() {
        given:
        def base = new MachineRequirementOpts([provisioning: 'spot'])

        when:
        def result = HintHelper.overlayHints(base, [consumableResources: 'my-license'])
        then:
        result.provisioning == 'spot'
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
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.maxSpotAttempts': '3'])
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
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.diskThroughputMiBps': '500'])
        then:
        result.diskThroughputMiBps == 500
    }

    def 'should overlay diskIops hint'() {
        given:
        def base = new MachineRequirementOpts([:])

        when:
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.diskIops': '10000'])
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
        def base = new MachineRequirementOpts([provisioning: 'ondemand'])

        when:
        def result = HintHelper.overlayHints(base, [
            'seqera/machineRequirement.provisioning': 'spotFirst',
            'seqera/machineRequirement.maxSpotAttempts': '3',
            'seqera/machineRequirement.diskType': 'ebs/gp3'
        ])
        then:
        result.provisioning == 'spotFirst'
        result.maxSpotAttempts == 3
        result.diskType == 'ebs/gp3'
    }

    def 'should preserve base values not overridden by hints'() {
        given:
        def base = new MachineRequirementOpts([provisioning: 'spot', diskType: 'ebs/gp3', diskMountPath: '/data'])

        when:
        def result = HintHelper.overlayHints(base, ['seqera/machineRequirement.diskType': 'ebs/io1'])
        then:
        result.provisioning == 'spot'
        result.diskType == 'ebs/io1'
        result.diskMountPath == '/data'
    }

    def 'should derive known keys from MachineRequirementOpts declared fields'() {
        expect: 'KNOWN_KEYS covers every declared field of MachineRequirementOpts'
        HintHelper.KNOWN_KEYS.size() > 0
        for( final field : MachineRequirementOpts.declaredFields ) {
            if( field.synthetic || java.lang.reflect.Modifier.isStatic(field.modifiers) || field.name.startsWith('$') || field.name == 'metaClass' )
                continue
            assert HintHelper.KNOWN_KEYS.contains("machineRequirement.${field.name}".toString())
        }
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
            'seqera/machineRequirement.provisioning': 'spot'
        ])
        then:
        result.size() == 1
        result['machineRequirement.provisioning'] == 'spot'
    }

    def 'should handle null hints map'() {
        when:
        def result = HintHelper.extractSeqeraHints(null)
        then:
        result.isEmpty()
    }

    def 'should accept unprefixed known keys'() {
        when:
        def result = HintHelper.extractSeqeraHints([
            'machineRequirement.provisioning': 'spot',
            'machineRequirement.diskType': 'ebs/gp3',
        ])
        then:
        result['machineRequirement.provisioning'] == 'spot'
        result['machineRequirement.diskType'] == 'ebs/gp3'
    }

    def 'should give prefixed form precedence over unprefixed'() {
        when:
        def result = HintHelper.extractSeqeraHints([
            'machineRequirement.provisioning': 'ondemand',
            'seqera/machineRequirement.provisioning': 'spotFirst',
        ])
        then:
        result['machineRequirement.provisioning'] == 'spotFirst'
    }

    def 'should overlay unprefixed hint onto base opts'() {
        given:
        def base = new MachineRequirementOpts([provisioning: 'ondemand'])

        when:
        def result = HintHelper.overlayHints(base, ['machineRequirement.provisioning': 'spotFirst'])
        then:
        result.provisioning == 'spotFirst'
    }

    def 'should ignore unknown unprefixed keys'() {
        when:
        def result = HintHelper.extractSeqeraHints([
            consumableResources: 'license-a=1',
            somethingElse: 'x',
            'machineRequirement.provisioning': 'spot',
        ])
        then:
        result.size() == 1
        result['machineRequirement.provisioning'] == 'spot'
    }

}
