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

import nextflow.cloud.types.PriceModel
import software.amazon.awssdk.services.batch.BatchClient
import software.amazon.awssdk.services.ec2.model.Instance
import software.amazon.awssdk.services.ec2.model.InstanceType
import software.amazon.awssdk.services.ec2.model.InstanceLifecycleType
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

    def 'should map unknown generation 8 instance types to sdk sentinel'() {
        expect:
        InstanceType.fromValue('r8id.xlarge') == InstanceType.UNKNOWN_TO_SDK_VERSION
        InstanceType.fromValue('m8id.xlarge') == InstanceType.UNKNOWN_TO_SDK_VERSION
        InstanceType.fromValue('c8id.xlarge') == InstanceType.UNKNOWN_TO_SDK_VERSION
    }
}
