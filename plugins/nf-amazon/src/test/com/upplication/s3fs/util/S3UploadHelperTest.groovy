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

package com.upplication.s3fs.util


import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class S3UploadHelperTest extends Specification {

    @Shared final long _1_KiB = 1_024
    @Shared final long _1_MiB = _1_KiB **2
    @Shared final long _1_GiB = _1_KiB **3
    @Shared final long _1_TiB = _1_KiB **4
    
    @Shared final long _10_MiB = _1_MiB * 10
    @Shared final long _100_MiB = _1_MiB * 100

    @Unroll
    def 'should compute s3 file chunk size' () {

        expect:
        S3UploadHelper.computePartSize(FILE_SIZE, CHUNK_SIZE) == EXPECTED_CHUNK_SIZE
        and:
        def parts = FILE_SIZE / EXPECTED_CHUNK_SIZE
        parts <= S3UploadHelper.MAX_PARTS_COUNT
        parts > 0

        where:
        FILE_SIZE   | EXPECTED_CHUNK_SIZE   | CHUNK_SIZE
        _1_KiB      | _10_MiB               | _10_MiB
        _1_MiB      | _10_MiB               | _10_MiB
        _1_GiB      | _10_MiB               | _10_MiB
        _1_TiB      | 110 * _1_MiB          | _10_MiB
        5 * _1_TiB  | 530 * _1_MiB          | _10_MiB
        10 * _1_TiB | 1050 * _1_MiB         | _10_MiB
        and:
        _1_KiB      | _100_MiB              | _100_MiB
        _1_MiB      | _100_MiB              | _100_MiB
        _1_GiB      | _100_MiB              | _100_MiB
        _1_TiB      | 110 * _1_MiB          | _100_MiB
        5 * _1_TiB  | 530 * _1_MiB          | _100_MiB
        10 * _1_TiB | 1050 * _1_MiB         | _100_MiB

    }


    def 'should check s3 part size' () {
        when:
        S3UploadHelper.checkPartSize(S3UploadHelper.MIN_PART_SIZE)
        then:
        noExceptionThrown()

        when:
        S3UploadHelper.checkPartSize(S3UploadHelper.MIN_PART_SIZE+1)
        then:
        noExceptionThrown()

        when:
        S3UploadHelper.checkPartSize(S3UploadHelper.MAX_PART_SIZE-1)
        then:
        noExceptionThrown()

        when:
        S3UploadHelper.checkPartSize(S3UploadHelper.MAX_PART_SIZE)
        then:
        noExceptionThrown()

        when:
        S3UploadHelper.checkPartSize(S3UploadHelper.MAX_PART_SIZE+1)
        then:
        thrown(IllegalArgumentException)

        when:
        S3UploadHelper.checkPartSize(S3UploadHelper.MIN_PART_SIZE-1)
        then:
        thrown(IllegalArgumentException)
    }

    def 'should check part index' () {
        when:
        S3UploadHelper.checkPartIndex(1, 's3://foo', 1000, 100)
        then:
        noExceptionThrown()

        when:
        S3UploadHelper.checkPartIndex(S3UploadHelper.MAX_PARTS_COUNT, 's3://foo', 1000, 100)
        then:
        noExceptionThrown()

        when:
        S3UploadHelper.checkPartIndex(S3UploadHelper.MAX_PARTS_COUNT+1, 's3://foo', 1000, 100)
        then:
        def e1 = thrown(IllegalArgumentException)
        e1.message == "S3 multipart copy request exceed the number of max allowed parts -- offending value: 10001; file: 's3://foo'; size: 1000; part-size: 100"

        when:
        S3UploadHelper.checkPartIndex(0, 's3://foo', 1000, 100)
        then:
        def e2 = thrown(IllegalArgumentException)
        e2.message == "S3 multipart copy request index cannot less than 1 -- offending value: 0; file: 's3://foo'; size: 1000; part-size: 100"


    }
}
