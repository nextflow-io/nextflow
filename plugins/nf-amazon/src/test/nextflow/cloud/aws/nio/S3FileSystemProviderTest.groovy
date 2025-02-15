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
 *
 */

package nextflow.cloud.aws.nio

import nextflow.cloud.aws.config.AwsConfig
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class S3FileSystemProviderTest extends Specification {

    @Unroll
    def 'should get global region' () {
        given:
        def provider = Spy(S3FileSystemProvider)

        expect:
        provider.globalRegion(new AwsConfig(CONFIG)) == EXPECTED

        where:
        EXPECTED    | CONFIG
        'us-east-1' | [:]
        'us-east-1' | [region:'foo']
        'us-east-1' | [region:'foo', client:[endpoint: 'http://s3.us-east-2.amazonaws.com']]
        'foo'       | [region:'foo', client:[endpoint: 'http://bar.com']]        

    }

}
