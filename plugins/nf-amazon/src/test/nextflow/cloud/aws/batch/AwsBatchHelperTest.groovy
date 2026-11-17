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

package nextflow.cloud.aws.batch

import nextflow.cloud.aws.AwsClientFactory
import nextflow.cloud.types.CloudMachineInfo
import nextflow.cloud.types.PriceModel
import software.amazon.awssdk.services.batch.BatchClient
import software.amazon.awssdk.services.ec2.Ec2Client
import software.amazon.awssdk.services.ec2.model.DescribeInstancesRequest
import software.amazon.awssdk.services.ec2.model.DescribeInstancesResponse
import software.amazon.awssdk.services.ec2.model.Instance
import software.amazon.awssdk.services.ec2.model.InstanceType
import software.amazon.awssdk.services.ec2.model.InstanceLifecycleType
import software.amazon.awssdk.services.ec2.model.Placement
import software.amazon.awssdk.services.ec2.model.Reservation
import spock.lang.Specification
import spock.lang.Unroll

/**
 * Tests for AwsBatchHelper
 *
 * @author Rob Syme <rob.syme@seqera.io>
 */
class AwsBatchHelperTest extends Specification {

    @Unroll
    def 'should detect spot instance pricing model'() {
        given:
        def helper = new AwsBatchHelper(Mock(BatchClient), null)
        def instance = Instance.builder()
            .instanceLifecycle(LIFECYCLE)
            .build()

        when:
        def result = helper.getPrice(instance)

        then:
        result == EXPECTED

        where:
        LIFECYCLE                       | EXPECTED
        InstanceLifecycleType.SPOT      | PriceModel.spot
        InstanceLifecycleType.SCHEDULED | PriceModel.standard
        null                            | PriceModel.standard   // on-demand instances return null
    }

    def 'should preserve raw aws instance type values'() {
        given:
        def helper = new AwsBatchHelper(Mock(BatchClient), null)

        expect:
        helper.getInstanceType(INSTANCE) == TYPE

        where:
        TYPE            | _
        'm4.large'      | _
        'r8id.xlarge'   | _
        and:
        INSTANCE = Instance.builder().instanceType(InstanceType.fromValue(TYPE)).instanceType(TYPE).build()
    }

    def 'should get machine info by instance id'() {
        given:
        def instance = Instance.builder()
            .instanceId('i-0abc123')
            .instanceType('m5.large')
            .placement(Placement.builder().availabilityZone('us-east-1a').build())
            .instanceLifecycle(InstanceLifecycleType.SPOT)
            .build()
        def requests = []
        def ec2 = new Ec2Client() {
            String serviceName() { 'ec2' }
            void close() {}
            DescribeInstancesResponse describeInstances(DescribeInstancesRequest req) {
                requests << req
                DescribeInstancesResponse.builder()
                    .reservations(Reservation.builder().instances(instance).build())
                    .build()
            }
        }
        def factory = Mock(AwsClientFactory) { getEc2Client() >> ec2 }
        def helper = new AwsBatchHelper(Mock(BatchClient), factory)

        when:
        def result = helper.getInfoByInstanceId('i-0abc123')

        then:
        requests*.instanceIds() == [['i-0abc123']]
        result == new CloudMachineInfo('m5.large', 'us-east-1a', PriceModel.spot, 'i-0abc123')
    }
}
