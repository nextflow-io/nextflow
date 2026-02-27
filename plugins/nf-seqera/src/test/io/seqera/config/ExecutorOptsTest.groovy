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

import nextflow.util.Duration
import spock.lang.Specification

/**
 * Unit tests for SeqeraExecutorConfig
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ExecutorOptsTest extends Specification {

    def 'should throw error when endpoint is missing' () {
        when:
        new ExecutorOpts([:])

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('Missing Seqera endpoint')
    }

    def 'should create config with minimal settings' () {
        when:
        def config = new ExecutorOpts([endpoint: 'https://sched.example.com'])

        then:
        config.endpoint == 'https://sched.example.com'
        config.region == 'eu-central-1'  // default
        config.keyPairName == null
        config.batchFlushInterval == Duration.of('1 sec')
        config.machineRequirement != null
        config.machineRequirement.arch == null
        config.machineRequirement.provisioning == null
        !config.autoLabels
    }

    def 'should create config with custom region' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            region: 'us-west-2'
        ])

        then:
        config.endpoint == 'https://sched.example.com'
        config.region == 'us-west-2'
    }

    def 'should create config with custom batch flush interval' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            batchFlushInterval: '5 sec'
        ])

        then:
        config.batchFlushInterval == Duration.of('5 sec')
    }

    def 'should create config with machine requirement settings' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            machineRequirement: [
                arch: 'arm64',
                provisioning: 'spotFirst',
                maxSpotAttempts: 3,
                machineTypes: ['m6g', 'c6g']
            ]
        ])

        then:
        config.machineRequirement != null
        config.machineRequirement.arch == 'arm64'
        config.machineRequirement.provisioning == 'spotFirst'
        config.machineRequirement.maxSpotAttempts == 3
        config.machineRequirement.machineTypes == ['m6g', 'c6g']
    }

    def 'should create config with retry policy' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            retryPolicy: [maxAttempts: 5, delay: '2s']
        ])

        then:
        config.retryOpts().maxAttempts == 5
        config.retryOpts().delay == Duration.of('2s')
    }

    def 'should create config with all settings' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            region: 'eu-west-1',
            keyPairName: 'my-key',
            batchFlushInterval: '2 sec',
            machineRequirement: [
                arch: 'x86_64',
                provisioning: 'spot'
            ]
        ])

        then:
        config.endpoint == 'https://sched.example.com'
        config.region == 'eu-west-1'
        config.keyPairName == 'my-key'
        config.batchFlushInterval == Duration.of('2 sec')
        config.machineRequirement.arch == 'x86_64'
        config.machineRequirement.provisioning == 'spot'
    }

    def 'should create config with labels' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            labels: [
                project: 'genomics',
                team: 'research',
                costCenter: 'CC-1234'
            ]
        ])

        then:
        config.labels == [project: 'genomics', team: 'research', costCenter: 'CC-1234']
    }

    def 'should handle null labels' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com'
        ])

        then:
        config.labels == null
    }

    def 'should handle empty labels' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            labels: [:]
        ])

        then:
        config.labels == [:]
    }

    def 'should enable auto labels' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            autoLabels: true
        ])

        then:
        config.autoLabels
    }

    def 'should create config with prediction model' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            predictionModel: 'qr/v1'
        ])

        then:
        config.predictionModel == 'qr/v1'
    }

    def 'should default prediction model to null' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com'
        ])

        then:
        config.predictionModel == null
    }

    def 'should create config with taskEnvironment' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            taskEnvironment: [FOO: 'bar', BAZ: 'qux']
        ])

        then:
        config.taskEnvironment == [FOO: 'bar', BAZ: 'qux']
    }

    def 'should handle null taskEnvironment' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com'
        ])

        then:
        config.taskEnvironment == null
    }

    def 'should handle empty taskEnvironment' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            taskEnvironment: [:]
        ])

        then:
        config.taskEnvironment == [:]
    }

    def 'should reject invalid prediction model' () {
        when:
        new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            predictionModel: 'invalid'
        ])

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains("Invalid prediction model 'invalid'")
        e.message.contains('qr/v1')
    }

}
