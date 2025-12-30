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

import spock.lang.Specification
import spock.lang.Unroll

/**
 * Tests for {@link AwsEcsTaskHandler}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsEcsTaskHandlerTest extends Specification {

    def 'should sanitize family name'() {
        given:
        def handler = [:] as AwsEcsTaskHandler

        expect:
        handler.sanitizeFamilyName('ubuntu:latest') == 'ubuntu-latest'
        handler.sanitizeFamilyName('quay.io/org/image:1.0') == 'quayio-org-image-10'
        handler.sanitizeFamilyName('registry.example.com/path/to/image@sha256:abc123') == 'registryexamplecom-path-to-image-sha256-abc123'
        handler.sanitizeFamilyName('simple') == 'simple'
        handler.sanitizeFamilyName(null) == 'unknown'
        handler.sanitizeFamilyName('') == 'unknown'
    }

    def 'should limit family name length'() {
        given:
        def handler = [:] as AwsEcsTaskHandler
        def longName = 'a' * 250

        when:
        def result = handler.sanitizeFamilyName(longName)

        then:
        result.size() == 200
    }

    @Unroll
    def 'should detect spot interruption: #stoppedReason'() {
        given:
        def handler = [:] as AwsEcsTaskHandler

        expect:
        handler.isSpotInterruption(stoppedReason, stopCode) == expected

        where:
        stoppedReason                                   | stopCode              | expected
        null                                            | null                  | false
        'Task finished successfully'                    | null                  | false
        'spot capacity unavailable'                     | null                  | true
        'Spot instance interruption'                    | null                  | true
        'Task stopped due to SPOT termination'          | null                  | true
        'Insufficient capacity in region'               | null                  | true
        'Normal termination'                            | 'SpotInterruption'    | true
        'User requested termination'                    | null                  | false
    }
}
