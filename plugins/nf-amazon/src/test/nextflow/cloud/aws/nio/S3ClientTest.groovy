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
package nextflow.cloud.aws.nio

import software.amazon.awssdk.services.s3.model.ChecksumMode
import software.amazon.awssdk.services.s3.model.HeadObjectRequest
import software.amazon.awssdk.services.s3.model.HeadObjectResponse
import spock.lang.Specification

class S3ClientTest extends Specification {

    def 'getObjectMetadata enables ChecksumMode so x-amz-checksum-* headers are returned'() {
        given:
        def aws = Mock(software.amazon.awssdk.services.s3.S3Client)
        def s3  = new S3Client(aws)
        HeadObjectRequest captured = null

        when:
        s3.getObjectMetadata('bucket', 'key')

        then:
        1 * aws.headObject(_ as HeadObjectRequest) >> { HeadObjectRequest req ->
            captured = req
            HeadObjectResponse.builder().build()
        }
        captured.checksumMode() == ChecksumMode.ENABLED
        captured.bucket() == 'bucket'
        captured.key() == 'key'
    }
}
