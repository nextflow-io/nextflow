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

package nextflow.cloud.aws

import spock.lang.Ignore
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AmazonClientFactoryTest extends Specification {

    @Ignore
    def 'should create client' () {
        when:
        def opts = [
                accessKey: '<SOME KEY>',
                secretKey: '<SOME SECRET>',
                region: 'eu-west-1',
                assumeRoleArn: '<SOME ROLE>'
        ]
        def client = new AmazonClientFactory(opts)
        and:
        def s3 = client.getS3Client()
        def buckets = s3.listBuckets()
        println buckets.name
        then:
        true
    }

}
