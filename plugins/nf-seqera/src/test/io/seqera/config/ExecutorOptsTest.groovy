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
        config.region == null
        config.provider == null
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

    def 'should enable all auto labels when set to true' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            autoLabels: true
        ])

        then:
        config.autoLabels == ExecutorOpts.VALID_AUTO_LABELS
    }

    def 'should disable auto labels when set to false' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            autoLabels: false
        ])

        then:
        config.autoLabels.isEmpty()
    }

    def 'should accept auto labels as a list of short names' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            autoLabels: ['runName', 'projectName']
        ])

        then:
        config.autoLabels == ['runName', 'projectName'] as Set
    }

    def 'should trim whitespace in auto labels list entries' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            autoLabels: [' runName', 'projectName ']
        ])

        then:
        config.autoLabels == ['runName', 'projectName'] as Set
    }

    def 'should accept auto labels as a comma-separated string' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            autoLabels: 'runName,projectName,workflowId'
        ])

        then:
        config.autoLabels == ['runName', 'projectName', 'workflowId'] as Set
    }

    def 'should tolerate whitespace around comma-separated auto labels' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            autoLabels: 'runName, projectName ,workflowId'
        ])

        then:
        config.autoLabels == ['runName', 'projectName', 'workflowId'] as Set
    }

    def 'should treat empty auto labels list as disabled' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            autoLabels: []
        ])

        then:
        config.autoLabels.isEmpty()
    }

    def 'should treat empty auto labels string as disabled' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            autoLabels: ''
        ])

        then:
        config.autoLabels.isEmpty()
    }

    def 'should reject unknown auto labels name' () {
        when:
        new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            autoLabels: ['runName', 'foo']
        ])

        then:
        def err = thrown(IllegalArgumentException)
        err.message.contains("'seqera.executor.autoLabels'")
        err.message.contains('foo')
        err.message.contains('valid names')
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

    def 'should create config with computeEnvId' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            computeEnvId: 'ce-12345'
        ])

        then:
        config.computeEnvId == 'ce-12345'
    }

    def 'should default computeEnvId to null' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com'
        ])

        then:
        config.computeEnvId == null
    }

    def 'should create config with provider' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            provider: 'aws'
        ])

        then:
        config.provider == 'aws'
    }

    def 'should default provider to null' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com'
        ])

        then:
        config.provider == null
    }

    def 'should create config with provider and region' () {
        when:
        def config = new ExecutorOpts([
            endpoint: 'https://sched.example.com',
            provider: 'aws',
            region: 'us-west-2'
        ])

        then:
        config.provider == 'aws'
        config.region == 'us-west-2'
    }


}
