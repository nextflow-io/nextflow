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

package com.upplication.s3fs.ng

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Paths
import java.time.temporal.ChronoUnit

import com.amazonaws.services.s3.AmazonS3
import com.amazonaws.services.s3.model.ObjectMetadata
import com.upplication.s3fs.S3FileSystem
import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.function.ContextualSupplier
import groovy.util.logging.Slf4j
import spock.lang.Ignore
import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Shared
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@IgnoreIf({System.getenv('NXF_SMOKE')})
@Requires({System.getenv('AWS_S3FS_ACCESS_KEY') && System.getenv('AWS_S3FS_SECRET_KEY')})
@Slf4j
class S3ParallelDownloadTest extends Specification {

    @Shared
    AmazonS3 s3Client

    def setupSpec() {
        def accessKey = System.getenv('AWS_S3FS_ACCESS_KEY') 
        def secretKey = System.getenv('AWS_S3FS_SECRET_KEY')
//        def region = System.getenv('AWS_REGION') ?: 'eu-west-1'
//        log.debug "Creating AWS S3 client: region=$region; accessKey=${accessKey?.substring(0,5)}.. - secretKey=${secretKey?.substring(0,5)}.. -  "
//        final creds = new AWSCredentials() {
//            String getAWSAccessKeyId() { accessKey }
//            String getAWSSecretKey() { secretKey }
//        }
//
//        storageClient = AmazonS3ClientBuilder
//                .standard()
//                .withRegion(region)
//                .withCredentials(new AWSStaticCredentialsProvider(creds))
//                .build()
        def fs = (S3FileSystem) FileSystems.newFileSystem(URI.create("s3:///"), [access_key: accessKey, secret_key: secretKey])
        s3Client = fs.client.getClient()
    }

    @Ignore
    def 'should download small file' () {
        given:
        def downloader = new S3ParallelDownload(s3Client)
        
        when:
        def stream = downloader.download('nextflow-ci','hello.txt')
        then:
        stream.text == 'Hello world\n'
    }

    @Ignore
    def 'should download 100 mb file'  () {
        given:
        def downloader = new S3ParallelDownload(s3Client)
        def target = Paths.get('file-100MB.data-copy.data')
        and:
        Files.deleteIfExists(target)
        when:
        def stream = downloader.download('nextflow-ci','file-100MB.data')
        then:
        Files.copy(stream, target)

        cleanup:
        stream?.close()
    }

    @Ignore
    def 'should download 10 gb file'  () {
        given:
        def downloader = new S3ParallelDownload(s3Client)
        def target = Paths.get('real.fastq.gz')

        when:
        Files.deleteIfExists(target)
        and:
        def stream = downloader.download('nextflow-ci','petagene/example_data/real.fastq.gz')
        then:
        Files.copy(stream, target)

        cleanup:
        stream?.close()
    }

    def 'should create part single' () {
        given:
        def FILE_LEN = 1
        def CHUNK_SIZE = 1000
        and:
        def client = Mock(AmazonS3)
        def download = new S3ParallelDownload(client, new DownloadOpts(download_chunk_size: String.valueOf(CHUNK_SIZE)))
        def META = Mock(ObjectMetadata) {getContentLength() >> FILE_LEN }

        when:
        def result = download.prepareGetPartRequests('foo','bar').iterator()
        then:
        1 * client.getObjectMetadata('foo','bar') >> META
        and:
        with(result.next()) {
            getBucketName() == 'foo'
            getKey() == 'bar'
            getRange() == [0,0]
        }
        and:
        !result.hasNext()
    }


    def 'should create part requests' () {
        given:
        def FILE_LEN = 3_000
        def CHUNK_SIZE = 1000
        and:
        def client = Mock(AmazonS3)
        def download = new S3ParallelDownload(client, new DownloadOpts(download_chunk_size: String.valueOf(CHUNK_SIZE)))
        def META = Mock(ObjectMetadata) {getContentLength() >> FILE_LEN }

        when:
        def result = download.prepareGetPartRequests('foo','bar').iterator()
        then:
        1 * client.getObjectMetadata('foo','bar') >> META
        and:
        with(result.next()) {
            getBucketName() == 'foo'
            getKey() == 'bar'
            getRange() == [0,999]
        }
        and:
        with(result.next()) {
            getBucketName() == 'foo'
            getKey() == 'bar'
            getRange() == [1000,1999]
        }
        and:
        with(result.next()) {
            getBucketName() == 'foo'
            getKey() == 'bar'
            getRange() == [2000,2999]
        }
        and:
        !result.hasNext()
    }

    def 'should create long requests' () {
        given:
        def FILE_LEN = 6_000_000_000
        def CHUNK_SIZE = 2_000_000_000
        and:
        def client = Mock(AmazonS3)
        def download = new S3ParallelDownload(client, new DownloadOpts(download_chunk_size: String.valueOf(CHUNK_SIZE), download_buffer_max_size: String.valueOf(CHUNK_SIZE)))
        def META = Mock(ObjectMetadata) {getContentLength() >> FILE_LEN }

        when:
        def result = download.prepareGetPartRequests('foo','bar')
        then:
        1 * client.getObjectMetadata('foo','bar') >> META
        and:
        with(result[0]) {
            getBucketName() == 'foo'
            getKey() == 'bar'
            getRange() == [0,1_999_999_999]
        }
        and:
        with(result[1]) {
            getBucketName() == 'foo'
            getKey() == 'bar'
            getRange() == [2_000_000_000,3_999_999_999]
        }
        and:
        with(result[2]) {
            getBucketName() == 'foo'
            getKey() == 'bar'
            getRange() == [4_000_000_000,5_999_999_999]
        }
        and:
        result.size()==3
    }

    @Ignore
    def 'test failsafe' () {
        given:
        RetryPolicy<Object> retryPolicy = RetryPolicy.builder()
                .handle(RuntimeException.class)
//                .withDelay(Duration.ofSeconds(1))
//                .withMaxDuration(Duration.of(60, ChronoUnit.SECONDS))
                .withBackoff(1, 30, ChronoUnit.SECONDS)
                .withMaxRetries(10)
                .onFailedAttempt(e -> log.error("Connection attempt failed - cause: ${e.getLastFailure()}"))
                .onRetry(e -> log.warn("Failure #{}. Retrying.", e.getAttemptCount()))
                .build();

        when:
        def work = { dev.failsafe.ExecutionContext it ->
            log.debug "try num ${it.getAttemptCount()}"
            throw new RuntimeException("Break ${it.getAttemptCount()}")
        } as ContextualSupplier
        def result = Failsafe.with(retryPolicy).get( work )
        then:
        result == 'Hello'
    }
}
