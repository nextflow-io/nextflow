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

package nextflow.cloud.google.batch.client

import nextflow.cloud.CloudTransferOptions
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BatchConfigTest extends Specification {

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
        !config.logsPath
        and:
        config.maxParallelTransfers == CloudTransferOptions.MAX_TRANSFER
        config.maxTransferAttempts == CloudTransferOptions.MAX_TRANSFER_ATTEMPTS
        config.delayBetweenAttempts == CloudTransferOptions.DEFAULT_DELAY_BETWEEN_ATTEMPTS
        !config.gcloudCli
        !config.gsutilCli
        and:
        !config.stageInCopyTransport
        !config.stageOutCopyTransport
        !config.usesGoogleBatchStaging()
    }

    def 'should reject invalid copy transport' () {
        when:
        new BatchConfig([stageInCopyTransport: 'ftp'])
        then:
        thrown(IllegalArgumentException)
    }

    def 'should detect google batch staging when cli transport set' () {
        expect:
        new BatchConfig([stageOutCopyTransport: 'gcloud']).usesGoogleBatchStaging()
        new BatchConfig([stageInCopyTransport: 'gsutil']).usesGoogleBatchStaging()
        !new BatchConfig([stageInCopyTransport: 'posix']).usesGoogleBatchStaging()
    }

    def 'should create batch config with custom settings' () {
        given:
        def opts = [
            spot: true,
            maxSpotAttempts: 8,
            autoRetryExitCodes: [50001, 50003, 50005],
            retryPolicy: [maxAttempts: 10],
            bootDiskImage: 'batch-foo',
            bootDiskSize: '100GB',
            logsPath: 'gs://my-logs-bucket/logs',
            installOpsAgent: true,
            stageInCopyTransport: 'gcloud',
            stageOutCopyTransport: 'posix',
            gcloudCli: '/opt/google/gcloud',
            maxParallelTransfers: 8,
            maxTransferAttempts: 3,
            delayBetweenAttempts: '5s'
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
        and:
        config.logsPath == 'gs://my-logs-bucket/logs'
        and:
        config.installOpsAgent == true
        and:
        config.stageInCopyTransport == 'gcloud'
        config.stageOutCopyTransport == 'posix'
        config.gcloudCli == '/opt/google/gcloud'
        config.maxParallelTransfers == 8
        config.maxTransferAttempts == 3
        config.delayBetweenAttempts == Duration.of('5s')
    }

}
