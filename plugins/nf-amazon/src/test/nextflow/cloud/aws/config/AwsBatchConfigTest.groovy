/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.cloud.aws.config

import java.nio.file.Paths

import nextflow.util.Duration
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsBatchConfigTest extends Specification {

    def 'should create default config' () {
        when:
        def batch = new AwsBatchConfig([:])
        then:
        batch.maxParallelTransfers == AwsBatchConfig.MAX_TRANSFER
        batch.maxTransferAttempts == AwsBatchConfig.DEFAULT_AWS_MAX_ATTEMPTS
        batch.delayBetweenAttempts == AwsBatchConfig.DEFAULT_DELAY_BETWEEN_ATTEMPTS
        batch.maxSpotAttempts == AwsBatchConfig.DEFAULT_MAX_SPOT_ATTEMPTS
        batch.retryMode == 'standard'
        and:
        !batch.cliPath
        !batch.volumes
        !batch.jobRole
        !batch.executionRole
        !batch.logsGroup
        !batch.shareIdentifier
        !batch.isFargateMode()
        !batch.s5cmdPath
        batch.schedulingPriority == 0
    }

    def 'should create config with options' () {
        given:
        def OPTS = [
                cliPath: '/some/bin/aws',
                maxParallelTransfers:1,
                maxTransferAttempts:2,
                delayBetweenAttempts: '3s',
                maxSpotAttempts: 4,
                volumes: '/some/path:/mnt/path,/other/path',
                jobRole: 'xyz',
                executionRole: 'some:exec:role',
                logsGroup: 'group-name-123',
                retryMode: 'legacy',
                shareIdentifier: 'id-x1',
                schedulingPriority: 100,
        ]

        when:
        def batch = new AwsBatchConfig(OPTS)
        then:
        batch.cliPath == '/some/bin/aws'
        batch.maxParallelTransfers == 1
        batch.maxTransferAttempts == 2
        batch.delayBetweenAttempts == Duration.of('3sec')
        batch.maxSpotAttempts == 4
        batch.volumes == ['/some/path:/mnt/path', '/other/path']
        batch.jobRole == 'xyz'
        batch.executionRole == 'some:exec:role'
        batch.logsGroup == 'group-name-123'
        batch.retryMode == 'legacy'
        batch.shareIdentifier == 'id-x1'
        batch.schedulingPriority == 100
        !batch.isFargateMode()
    }

    def 'should parse volumes list' () {

        given:
        def executor = Spy(AwsBatchConfig)

        expect:
        executor.makeVols(OBJ) == EXPECTED

        where:
        OBJ             | EXPECTED
        null            | []
        'foo'           | ['foo']
        'foo, bar'      | ['foo','bar']
        '/foo/,/bar///' | ['/foo','/bar']
        ['/this','/that'] | ['/this','/that']
        ['/foo/bar/']   | ['/foo/bar']

    }

    def 'should add a volume' () {
        given:
        def opts = new AwsBatchConfig()

        when:
        opts.addVolume(Paths.get('/some/dir'))
        then:
        opts.volumes == ['/some/dir']

        when:
        opts.addVolume(Paths.get('/other/dir'))
        opts.addVolume(Paths.get('/other/dir'))
        then:
        opts.volumes == ['/some/dir', '/other/dir']
    }

    def 'should parse cli path' () {
        given:
        def opts = new AwsBatchConfig(OPTS)

        expect:
        opts.cliPath == S3_CLI_PATH
        opts.s5cmdPath == S5_CLI_PATH
        opts.isFargateMode() == FARGATE
        
        where:
        OPTS                                                        | S3_CLI_PATH       | S5_CLI_PATH       | FARGATE
        [:]                                                         | null              | null              | false
        [cliPath: "/opt/bin/aws"]                                   | '/opt/bin/aws'    | null              | false
        [cliPath: "/s5cmd"]                                         | null              | null              | false
        [cliPath: "/opt/s5cmd --foo"]                               | null              | null              | false
        and:
        [platformType: 'fargate', cliPath: "/opt/bin/aws"]          | null              | 's5cmd'           | true
        [platformType: 'fargate', cliPath: "/opt/s5cmd"]            | null              | '/opt/s5cmd'      | true
        [platformType: 'fargate', cliPath: "/opt/s5cmd --foo"]      | null              | '/opt/s5cmd --foo'| true
    }

}
