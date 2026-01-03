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

import java.nio.file.Paths

import nextflow.Session
import nextflow.cloud.aws.nio.S3Path
import nextflow.exception.AbortOperationException
import nextflow.processor.TaskRun
import software.amazon.awssdk.services.ec2.Ec2Client
import software.amazon.awssdk.services.ec2.model.DescribeSecurityGroupsResponse
import software.amazon.awssdk.services.ec2.model.DescribeSubnetsResponse
import software.amazon.awssdk.services.ec2.model.DescribeVpcsResponse
import software.amazon.awssdk.services.ec2.model.SecurityGroup
import software.amazon.awssdk.services.ec2.model.Subnet
import software.amazon.awssdk.services.ec2.model.Vpc
import software.amazon.awssdk.services.ecs.EcsClient
import spock.lang.Specification

/**
 * Tests for {@link AwsEcsExecutor}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsEcsExecutorTest extends Specification {

    def 'should validate fusion is enabled'() {
        given:
        def session = Mock(Session) {
            getConfig() >> [aws: [ecs: [cluster: 'test-cluster']]]
        }
        def executor = Spy(AwsEcsExecutor)
        executor.session = session

        when:
        executor.isFusionEnabled() >> false
        executor.validateFusion()

        then:
        1 * session.abort()
        thrown(AbortOperationException)
    }

    def 'should pass fusion validation when enabled'() {
        given:
        def session = Mock(Session) {
            getConfig() >> [aws: [ecs: [cluster: 'test-cluster']]]
        }
        def executor = Spy(AwsEcsExecutor)
        executor.session = session

        when:
        executor.isFusionEnabled() >> true
        executor.validateFusion()

        then:
        noExceptionThrown()
    }

    def 'should validate work directory is S3'() {
        given:
        def session = Mock(Session) {
            getConfig() >> [aws: [ecs: [cluster: 'test-cluster']]]
        }
        def executor = Spy(AwsEcsExecutor) {
            getWorkDir() >> Paths.get('/local/work/dir')
        }
        executor.session = session

        when:
        executor.validateWorkDir()

        then:
        1 * session.abort()
        thrown(AbortOperationException)
    }

    def 'should pass work dir validation for S3 path'() {
        given:
        def s3Path = Mock(S3Path)
        def session = Mock(Session) {
            getConfig() >> [aws: [ecs: [cluster: 'test-cluster']]]
        }
        def executor = Spy(AwsEcsExecutor) {
            getWorkDir() >> s3Path
        }
        executor.session = session

        when:
        executor.validateWorkDir()

        then:
        noExceptionThrown()
    }

    def 'should use configured subnets'() {
        given:
        def awsOptions = Mock(AwsEcsOptions) {
            getSubnets() >> ['subnet-123', 'subnet-456']
            getSecurityGroups() >> null
        }
        def ec2Client = Mock(Ec2Client)
        def executor = new AwsEcsExecutor()
        executor.@awsOptions = awsOptions
        executor.@ec2Client = ec2Client

        when:
        executor.discoverVpcConfiguration()

        then:
        executor.resolvedSubnets == ['subnet-123', 'subnet-456']
        // Security groups should be auto-discovered
        1 * ec2Client.describeVpcs(_) >> DescribeVpcsResponse.builder()
            .vpcs(Vpc.builder().vpcId('vpc-default').build())
            .build()
        1 * ec2Client.describeSecurityGroups(_) >> DescribeSecurityGroupsResponse.builder()
            .securityGroups(SecurityGroup.builder().groupId('sg-default').build())
            .build()
    }

    def 'should use configured security groups'() {
        given:
        def awsOptions = Mock(AwsEcsOptions) {
            getSubnets() >> null
            getSecurityGroups() >> ['sg-123']
        }
        def ec2Client = Mock(Ec2Client)
        def executor = new AwsEcsExecutor()
        executor.@awsOptions = awsOptions
        executor.@ec2Client = ec2Client

        when:
        executor.discoverVpcConfiguration()

        then:
        executor.resolvedSecurityGroups == ['sg-123']
        // Subnets should be auto-discovered (one describeVpcs call for subnet discovery)
        1 * ec2Client.describeVpcs(_) >> DescribeVpcsResponse.builder()
            .vpcs(Vpc.builder().vpcId('vpc-default').build())
            .build()
        1 * ec2Client.describeSubnets(_) >> DescribeSubnetsResponse.builder()
            .subnets(
                Subnet.builder().subnetId('subnet-a').build(),
                Subnet.builder().subnetId('subnet-b').build()
            )
            .build()
    }

    def 'should auto-discover VPC configuration'() {
        given:
        def awsOptions = Mock(AwsEcsOptions) {
            getSubnets() >> null
            getSecurityGroups() >> null
        }
        def ec2Client = Mock(Ec2Client)
        def executor = new AwsEcsExecutor()
        executor.@awsOptions = awsOptions
        executor.@ec2Client = ec2Client

        when:
        executor.discoverVpcConfiguration()

        then:
        // Called twice: once for subnets, once for security groups
        2 * ec2Client.describeVpcs(_) >> DescribeVpcsResponse.builder()
            .vpcs(Vpc.builder().vpcId('vpc-default').build())
            .build()
        1 * ec2Client.describeSubnets(_) >> DescribeSubnetsResponse.builder()
            .subnets(
                Subnet.builder().subnetId('subnet-a').build(),
                Subnet.builder().subnetId('subnet-b').build()
            )
            .build()
        1 * ec2Client.describeSecurityGroups(_) >> DescribeSecurityGroupsResponse.builder()
            .securityGroups(SecurityGroup.builder().groupId('sg-default').build())
            .build()

        and:
        executor.resolvedSubnets == ['subnet-a', 'subnet-b']
        executor.resolvedSecurityGroups == ['sg-default']
    }

    def 'should fail when no default VPC exists for subnet discovery'() {
        given:
        def awsOptions = Mock(AwsEcsOptions) {
            getSubnets() >> null
            getSecurityGroups() >> ['sg-123']
        }
        def ec2Client = Mock(Ec2Client)
        def executor = new AwsEcsExecutor()
        executor.@awsOptions = awsOptions
        executor.@ec2Client = ec2Client

        when:
        executor.discoverVpcConfiguration()

        then:
        1 * ec2Client.describeVpcs(_) >> DescribeVpcsResponse.builder()
            .vpcs([])
            .build()

        thrown(AbortOperationException)
    }

    def 'should fail when no subnets found in default VPC'() {
        given:
        def awsOptions = Mock(AwsEcsOptions) {
            getSubnets() >> null
            getSecurityGroups() >> ['sg-123']
        }
        def ec2Client = Mock(Ec2Client)
        def executor = new AwsEcsExecutor()
        executor.@awsOptions = awsOptions
        executor.@ec2Client = ec2Client

        when:
        executor.discoverVpcConfiguration()

        then:
        1 * ec2Client.describeVpcs(_) >> DescribeVpcsResponse.builder()
            .vpcs(Vpc.builder().vpcId('vpc-default').build())
            .build()
        1 * ec2Client.describeSubnets(_) >> DescribeSubnetsResponse.builder()
            .subnets([])
            .build()

        thrown(AbortOperationException)
    }

    def 'should create task handler'() {
        given:
        def task = Mock(TaskRun) {
            getWorkDir() >> Paths.get('/work/dir')
        }
        def executor = new AwsEcsExecutor()
        executor.@awsOptions = Mock(AwsEcsOptions)
        executor.@ecsClient = Mock(EcsClient)

        when:
        def handler = executor.createTaskHandler(task)

        then:
        handler instanceof AwsEcsTaskHandler
    }

    def 'should return container native true'() {
        given:
        def executor = new AwsEcsExecutor()

        expect:
        executor.isContainerNative() == true
    }

    def 'should return secret native true'() {
        given:
        def executor = new AwsEcsExecutor()

        expect:
        executor.isSecretNative() == true
    }

    def 'should return docker container engine'() {
        given:
        def executor = new AwsEcsExecutor()

        expect:
        executor.containerConfigEngine() == 'docker'
    }

    def 'should warn about PATH env in config'() {
        given:
        def session = Mock(Session) {
            getConfig() >> [env: [PATH: '/some/path']]
        }
        def executor = new AwsEcsExecutor()
        executor.session = session

        when:
        executor.validatePathDir()

        then:
        // Should complete without error (just logs a warning)
        noExceptionThrown()
    }
}
