/*
 * Copyright 2013-2024, Seqera Labs
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
            .instanceLifecycle(lifecycle)
            .build()

        when:
        def result = helper.getPrice(instance)

        then:
        result == expected

        where:
        lifecycle                       | expected
        InstanceLifecycleType.SPOT      | PriceModel.spot
        InstanceLifecycleType.SCHEDULED | PriceModel.standard
        null                            | PriceModel.standard   // on-demand instances return null
    }
}
