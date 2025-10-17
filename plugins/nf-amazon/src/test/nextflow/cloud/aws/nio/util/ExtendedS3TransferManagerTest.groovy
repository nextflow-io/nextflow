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

package nextflow.cloud.aws.nio.util

import software.amazon.awssdk.transfer.s3.S3TransferManager
import spock.lang.Specification
import spock.lang.Unroll


/**
 * Test for ExtendedS3TransferManager
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ExtendedS3TransferManagerTest extends Specification {

    def 'should initialize with default values'() {
        given:
        def mockTransferManager = Mock(S3TransferManager)
        def props = new Properties()

        when:
        def extendedManager = new ExtendedS3TransferManager(mockTransferManager, props)

        then:
        extendedManager.partSize == 8 * 1024 * 1024 // 8 MB
        extendedManager.downloadPermits == 50 // 400MB / 8MB
    }

    def 'should initialize with custom properties'() {
        given:
        def mockTransferManager = Mock(S3TransferManager)
        def props = new Properties()
        props.setProperty('max_download_heap_memory', '200000000') // 200 MB
        props.setProperty('minimum_part_size', '16777216')    // 16 MB

        when:
        def extendedManager = new ExtendedS3TransferManager(mockTransferManager, props)

        then:
        extendedManager.partSize == 16 * 1024 * 1024 // 16 MB
        extendedManager.downloadPermits == 11 // 200MB / 16MB (floor) = 11.92... -> 11
    }

    @Unroll
    def 'should estimate parts correctly'() {
        given:
        def mockTransferManager = Mock(S3TransferManager)
        def props = new Properties()
        props.setProperty('minimum_part_size', partSizeStr)
        def extendedManager = new ExtendedS3TransferManager(mockTransferManager, props)

        expect:
        extendedManager.estimateParts(fileSize) == expectedParts

        where:
        fileSize        | partSizeStr    | expectedParts
        1024           | '8388608'      | 1           // 1KB file, 8MB parts = 1 part
        8388608        | '8388608'      | 1           // 8MB file, 8MB parts = 1 part
        16777216       | '8388608'      | 2           // 16MB file, 8MB parts = 2 parts
        100000000      | '8388608'      | 10          // ~95MB file, 8MB parts = 10 parts (capped at DEFAULT_INIT_BUFFER_PARTS)
        500000000      | '8388608'      | 10          // ~476MB file, 8MB parts = 10 parts (capped at DEFAULT_INIT_BUFFER_PARTS)
        1048576        | '1048576'      | 1           // 1MB file, 1MB parts = 1 part
        10485760       | '1048576'      | 10          // 10MB file, 1MB parts = 10 parts (capped at DEFAULT_INIT_BUFFER_PARTS)
    }


    def 'should calculate downloadPermits correctly'() {
        given:
        def mockTransferManager = Mock(S3TransferManager)
        def props = new Properties()
        props.setProperty('max_download_heap_memory', maxBuffer)
        props.setProperty('minimum_part_size', partSize)
        
        when:
        def extendedManager = new ExtendedS3TransferManager(mockTransferManager, props)

        then:
        extendedManager.downloadPermits == expectedMaxParts

        where:
        maxBuffer     | partSize     | expectedMaxParts
        '419430400'   | '8388608'    | 50              // 400MB / 8MB
        '104857600'   | '8388608'    | 12              // 100MB / 8MB
        '838860800'   | '8388608'    | 100             // 800MB / 8MB
        '419430400'   | '16777216'   | 25              // 400MB / 16MB
    }

    def 'should handle zero or negative file sizes in estimateParts'() {
        given:
        def mockTransferManager = Mock(S3TransferManager)
        def props = new Properties()
        def extendedManager = new ExtendedS3TransferManager(mockTransferManager, props)

        expect:
        extendedManager.estimateParts(0) == 1
        extendedManager.estimateParts(-100) == 1
    }

    def 'should handle large file sizes in estimateParts'() {
        given:
        def mockTransferManager = Mock(S3TransferManager)
        def props = new Properties()
        def extendedManager = new ExtendedS3TransferManager(mockTransferManager, props)

        when:
        def parts = extendedManager.estimateParts(Long.MAX_VALUE)

        then:
        parts == 10 // Should be capped at DEFAULT_INIT_BUFFER_PARTS
    }
}