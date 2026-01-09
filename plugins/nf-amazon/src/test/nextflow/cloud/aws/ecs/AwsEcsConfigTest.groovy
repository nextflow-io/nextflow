/*
 * Copyright 2020-2024, Seqera Labs
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

package nextflow.cloud.aws.ecs

import nextflow.cloud.aws.config.AwsEcsConfig
import spock.lang.Specification

/**
 * Tests for {@link AwsEcsConfig}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsEcsConfigTest extends Specification {

    def 'should create default config'() {
        when:
        def ecs = new AwsEcsConfig([:])

        then:
        ecs.cluster == null
        ecs.executionRole == null
        ecs.taskRole == null
        ecs.subnets == null
        ecs.securityGroups == null
        ecs.logsGroup == AwsEcsConfig.DEFAULT_LOGS_GROUP
        ecs.maxSpotAttempts == AwsEcsConfig.DEFAULT_MAX_SPOT_ATTEMPTS
        ecs.assignPublicIp == true
    }

    def 'should create config with options'() {
        given:
        def OPTS = [
            cluster: 'my-cluster',
            executionRole: 'arn:aws:iam::123456789:role/ecsTaskExecutionRole',
            taskRole: 'arn:aws:iam::123456789:role/ecsTaskRole',
            subnets: ['subnet-abc123', 'subnet-def456'],
            securityGroups: ['sg-xyz789'],
            logsGroup: '/custom/logs/group',
            maxSpotAttempts: 10,
            assignPublicIp: false
        ]

        when:
        def ecs = new AwsEcsConfig(OPTS)

        then:
        ecs.cluster == 'my-cluster'
        ecs.executionRole == 'arn:aws:iam::123456789:role/ecsTaskExecutionRole'
        ecs.taskRole == 'arn:aws:iam::123456789:role/ecsTaskRole'
        ecs.subnets == ['subnet-abc123', 'subnet-def456']
        ecs.securityGroups == ['sg-xyz789']
        ecs.logsGroup == '/custom/logs/group'
        ecs.maxSpotAttempts == 10
        ecs.assignPublicIp == false
    }

    def 'should parse string lists'() {
        given:
        def ecs = new AwsEcsConfig([:])

        expect:
        ecs.parseStringList(OBJ) == EXPECTED

        where:
        OBJ                         | EXPECTED
        null                        | null
        'foo'                       | ['foo']
        'foo, bar'                  | ['foo', 'bar']
        'subnet-1,subnet-2,subnet-3'| ['subnet-1', 'subnet-2', 'subnet-3']
        ['subnet-1', 'subnet-2']    | ['subnet-1', 'subnet-2']
    }

    def 'should use default values when not specified'() {
        when:
        def ecs = new AwsEcsConfig([
            cluster: 'test-cluster',
            executionRole: 'arn:aws:iam::123:role/test'
        ])

        then:
        ecs.cluster == 'test-cluster'
        ecs.executionRole == 'arn:aws:iam::123:role/test'
        ecs.logsGroup == '/aws/ecs/nextflow'
        ecs.maxSpotAttempts == 5
        ecs.assignPublicIp == true
    }

    def 'should handle maxSpotAttempts zero value'() {
        when:
        def ecs = new AwsEcsConfig([maxSpotAttempts: 0])

        then:
        ecs.maxSpotAttempts == 0
    }

    def 'should handle assignPublicIp false value'() {
        when:
        def ecs = new AwsEcsConfig([assignPublicIp: false])

        then:
        ecs.assignPublicIp == false
    }
}
