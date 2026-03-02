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

package io.seqera.config

import spock.lang.Specification

/**
 * Unit tests for MachineRequirementOpts
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class MachineRequirementOptsTest extends Specification {

    def 'should create with empty config' () {
        when:
        def opts = new MachineRequirementOpts([:])

        then:
        opts.arch == null
        opts.provisioning == null
        opts.maxSpotAttempts == null
        opts.machineTypes == null
    }

    def 'should create with all settings' () {
        when:
        def opts = new MachineRequirementOpts([
            arch: 'arm64',
            provisioning: 'spotFirst',
            maxSpotAttempts: 3,
            machineTypes: ['m5', 'c5', 'r5']
        ])

        then:
        opts.arch == 'arm64'
        opts.provisioning == 'spotFirst'
        opts.maxSpotAttempts == 3
        opts.machineTypes == ['m5', 'c5', 'r5']
    }

    def 'should create with partial settings' () {
        when:
        def opts = new MachineRequirementOpts([arch: 'x86_64', provisioning: 'spot'])

        then:
        opts.arch == 'x86_64'
        opts.provisioning == 'spot'
        opts.maxSpotAttempts == null
        opts.machineTypes == null
    }

}
