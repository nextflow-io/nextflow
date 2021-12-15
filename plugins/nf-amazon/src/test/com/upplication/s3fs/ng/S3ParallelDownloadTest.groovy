/*
 * Copyright 2020-2021, Seqera Labs
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

import com.amazonaws.services.s3.AmazonS3
import com.upplication.s3fs.S3FileSystem
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


    def 'should create part requests' () {
        given:
        def download = new S3ParallelDownload(Mock(AmazonS3), new DownloadOpts(download_chunk_size: '1000'))

        when:
        def result = download.prepareGetPartRequests('foo','bar', 3_000)
        then:
        result.size()==3

    }

    def 'validate priority index' () {
        expect:
        S3ParallelDownload.seqOrder(0,0) == 0
        S3ParallelDownload.seqOrder(0,1) == 1
        S3ParallelDownload.seqOrder(0,2) == 2
        and:
        S3ParallelDownload.seqOrder(1,0) == 65536
        S3ParallelDownload.seqOrder(1,1) == 65537
        S3ParallelDownload.seqOrder(1,2) == 65538
        and:
        S3ParallelDownload.seqOrder(2,0) == (65536 * 2)
        S3ParallelDownload.seqOrder(2,1) == (65536 * 2) +1
        S3ParallelDownload.seqOrder(2,2) == (65536 * 2) +2
    }

    double logistic(float x, double A, double K) {
        return A / ( 1 + Math.exp( -1 * K * x ) )
    }

    int delay( int current, int capacity ) {
        def x = current / capacity * 100
        Math.round(logistic(x-90, 100_000, 0.5))
    }

    def 'should compute buffer delay' () {
        when:
        120.times {
            println "f($it) => " + delay(it, 100)
        }
        then:
        noExceptionThrown()
    }
}
