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

package io.seqera.config

import nextflow.util.Duration
import spock.lang.Specification

/**
 * Unit tests for SeqeraConfig
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SeqeraConfigTest extends Specification {

    def 'should create config with no executor' () {
        when:
        def config = new SeqeraConfig([:])

        then:
        config.executor == null
    }

    def 'should create config with executor settings' () {
        when:
        def config = new SeqeraConfig([
            executor: [
                endpoint: 'https://sched.example.com',
                region: 'us-west-2'
            ]
        ])

        then:
        config.executor != null
        config.executor.endpoint == 'https://sched.example.com'
        config.executor.region == 'us-west-2'
    }

    def 'should create config with full executor settings' () {
        when:
        def config = new SeqeraConfig([
            executor: [
                endpoint: 'https://sched.example.com',
                region: 'eu-west-1',
                keyPairName: 'my-key',
                batchFlushInterval: '2 sec',
                machineRequirement: [
                    arch: 'arm64',
                    provisioning: 'spot'
                ],
                labels: [
                    project: 'genomics',
                    team: 'research'
                ]
            ]
        ])

        then:
        config.executor != null
        config.executor.endpoint == 'https://sched.example.com'
        config.executor.region == 'eu-west-1'
        config.executor.keyPairName == 'my-key'
        config.executor.batchFlushInterval == Duration.of('2 sec')
        config.executor.machineRequirement.arch == 'arm64'
        config.executor.machineRequirement.provisioning == 'spot'
        config.executor.labels == [project: 'genomics', team: 'research']
    }

    def 'should throw error when executor endpoint is missing' () {
        when:
        new SeqeraConfig([
            executor: [
                region: 'us-west-2'
            ]
        ])

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('Missing Seqera endpoint')
    }

}
