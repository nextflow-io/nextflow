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

package com.upplication.s3fs.experiment

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Paths

import com.amazonaws.services.s3.AmazonS3
import com.amazonaws.services.s3.transfer.TransferManagerBuilder
import com.upplication.s3fs.AwsS3BaseSpec
import com.upplication.s3fs.S3FileSystem
import spock.lang.Ignore
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Ignore
class FileTransferTest extends Specification implements AwsS3BaseSpec {

    @Shared
    AmazonS3 s3Client

    def setupSpec() {
        def accessKey = System.getenv('AWS_S3FS_ACCESS_KEY')
        def secretKey = System.getenv('AWS_S3FS_SECRET_KEY')
        def fs = (S3FileSystem) FileSystems.newFileSystem(URI.create("s3:///"), [access_key: accessKey, secret_key: secretKey])
        s3Client = fs.client.getClient()
    }

    def 'should upload big file' () {
        given:
        def bucket = createBucket()
        def uri = new URI("s3:///$bucket/big.file")
        def target = Paths.get(uri)
        log.debug "Should create S3 file=$uri"

        when:
        def big = Paths.get('/Users/pditommaso/Projects/nextflow/big.file')
        Files.copy(big, target)
        then:
        existsPath("$bucket/file-name.txt")

        cleanup:
        if( bucket ) deleteBucket(bucket)
    }


    def 'should transfer file' () {
        given:
        def manager = TransferManagerBuilder
                .standard()
                .withS3Client(s3Client)
                .build();
        def target = new File('big.file.txt')
        when:
        def download = manager.download('nextflow-ci', 'petagene/example_data/real.fastq.gz', target)
        and:
        download.waitForCompletion()
        then:
        noExceptionThrown()
    }

    def 'should download big file' () {
        given:
        def uri = new URI('s3:///nextflow-ci/petagene/example_data/real.fastq.gz')
        def target = Paths.get('real.fastq.gz')
        Files.deleteIfExists(target)
        when:
        println "Starting download"
        Files.copy(Paths.get(uri), target)
        then:
        noExceptionThrown()

        cleanup:
        Files.deleteIfExists(target)
    }

}
