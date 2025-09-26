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

package nextflow.cloud.google.batch.client

import nextflow.Session
import nextflow.util.MemoryUnit
import spock.lang.Requires
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BatchConfigTest extends Specification {

    @Requires({System.getenv('GOOGLE_APPLICATION_CREDENTIALS')})
    def 'should create batch config' () {
        when:
        def config = new BatchConfig([:])
        then:
        !config.getSpot()
        and:
        config.retryConfig.maxAttempts == 5
        config.maxSpotAttempts == 0
        config.autoRetryExitCodes == [50001]
        and:
        !config.bootDiskImage
        !config.bootDiskSize
    }

    @Requires({System.getenv('GOOGLE_APPLICATION_CREDENTIALS')})
    def 'should create batch config with custom settings' () {
        given:
        def opts = [
            spot: true,
            maxSpotAttempts: 8,
            autoRetryExitCodes: [50001, 50003, 50005],
            retryPolicy: [maxAttempts: 10],
            bootDiskImage: 'batch-foo',
            bootDiskSize: '100GB'
        ]

        when:
        def config = new BatchConfig(opts)
        then:
        config.getSpot()
        and:
        config.retryConfig.maxAttempts == 10
        config.maxSpotAttempts == 8
        config.autoRetryExitCodes == [50001, 50003, 50005]
        and:
        config.bootDiskImage == 'batch-foo'
        config.bootDiskSize == MemoryUnit.of('100GB')
    }

    @Requires({System.getenv('GOOGLE_APPLICATION_CREDENTIALS')})
    def 'should validate logs bucket config' () {
        when:
        def config = new BatchConfig([logsBucket: 'gs://my-logs-bucket/logs'])
        then:
        config.logsBucket == 'gs://my-logs-bucket/logs'

        when:
        config = new BatchConfig([:])
        then:
        config.logsBucket == null
    }

    @Requires({System.getenv('GOOGLE_APPLICATION_CREDENTIALS')})
    def 'should reject invalid logs bucket paths' () {
        when:
        new BatchConfig([logsBucket: 'invalid-bucket'])
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains("Logs bucket path must start with 'gs://'")

        when:
        new BatchConfig([logsBucket: 'gs://'])
        then:
        e = thrown(IllegalArgumentException)
        e.message.contains("Invalid logs bucket path")

        when:
        new BatchConfig([logsBucket: 's3://bucket'])
        then:
        e = thrown(IllegalArgumentException)
        e.message.contains("Logs bucket path must start with 'gs://'")
    }

    def 'should extract bucket name from GCS path' () {
        expect:
        BatchConfig.extractBucketName('gs://my-bucket') == 'my-bucket'
        BatchConfig.extractBucketName('gs://my-bucket/logs') == 'my-bucket'
        BatchConfig.extractBucketName('gs://my-bucket/path/to/logs') == 'my-bucket'
        BatchConfig.extractBucketName('gs://') == ''
        BatchConfig.extractBucketName('invalid-path') == null
        BatchConfig.extractBucketName(null) == null
    }

    def 'should convert GCS path to mount path' () {
        expect:
        BatchConfig.convertGcsPathToMountPath('gs://my-bucket') == '/mnt/disks/my-bucket'
        BatchConfig.convertGcsPathToMountPath('gs://my-bucket/logs') == '/mnt/disks/my-bucket/logs'
        BatchConfig.convertGcsPathToMountPath('gs://my-bucket/path/to/logs') == '/mnt/disks/my-bucket/path/to/logs'
        BatchConfig.convertGcsPathToMountPath('invalid-path') == 'invalid-path'
        BatchConfig.convertGcsPathToMountPath(null) == null
    }

}
