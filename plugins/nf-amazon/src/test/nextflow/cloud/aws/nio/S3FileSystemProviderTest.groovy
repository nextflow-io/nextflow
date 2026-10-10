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

import java.time.Instant

import software.amazon.awssdk.services.s3.model.ListMultipartUploadsRequest
import software.amazon.awssdk.services.s3.model.ListMultipartUploadsResponse
import software.amazon.awssdk.services.s3.model.ListPartsRequest
import software.amazon.awssdk.services.s3.model.ListPartsResponse
import software.amazon.awssdk.services.s3.model.MultipartUpload
import software.amazon.awssdk.services.s3.model.ObjectCannedACL
import software.amazon.awssdk.services.s3.model.Part
import software.amazon.awssdk.services.s3.model.S3Exception
import software.amazon.awssdk.services.s3.model.ServerSideEncryption
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class S3FileSystemProviderTest extends Specification {

    def 'should create filesystem from config'(){
        given:
        def config = [
            client: [
                anonymous: true,
                s3Acl: 'Private',
                connectionTimeout: 20000,
                endpoint: 'https://s3.eu-west-1.amazonaws.com',
                maxConcurrency: 10,
                maxNativeMemory: '500MB',
                minimumPartSize: '7MB',
                multipartThreshold: '32MB',
                maxConnections: 100,
                maxErrorRetry: 3,
                socketTimeout: 20000,
                requesterPays: true,
                s3PathStyleAccess: true,
                proxyHost: 'host.com',
                proxyPort: 80,
                proxyScheme: 'https',
                proxyUsername: 'user',
                proxyPassword: 'pass',
                storageEncryption: 'AES256',
                storageKmsKeyId: 'arn:key:id',
                uploadMaxThreads: 15,
                uploadChunkSize: '7MB',
                uploadMaxAttempts: 4,
                uploadRetrySleep: '200ms'
            ],
            accessKey: '123456abc',
            secretKey: '78910def',
            profile: 'test'
        ]
        def provider = new S3FileSystemProvider();
        when:
        def fs = provider.newFileSystem(new URI("s3:///bucket/key"), config) as S3FileSystem
        then:
        fs.getBucketName() == 'bucket'
        def client = fs.getClient()
        client.client != null
        client.cannedAcl == ObjectCannedACL.PRIVATE
        client.storageEncryption == ServerSideEncryption.AES256
        client.isRequesterPaysEnabled == true
        client.kmsKeyId == 'arn:key:id'
        client.factory.accessKey() == '123456abc'
        client.factory.secretKey() == '78910def'
        client.factory.profile() == 'test'
        client.factory.config.s3Config.anonymous == true
        client.factory.config.s3Config.endpoint == 'https://s3.eu-west-1.amazonaws.com'
        client.factory.config.s3Config.pathStyleAccess == true
        fs.properties().getProperty('proxy_host') == 'host.com'
        fs.properties().getProperty('proxy_port') == '80'
        fs.properties().getProperty('proxy_scheme') == 'https'
        fs.properties().getProperty('proxy_username') == 'user'
        fs.properties().getProperty('proxy_password') == 'pass'
        fs.properties().getProperty('socket_timeout') == '20000'
        fs.properties().getProperty('connection_timeout') == '20000'
        fs.properties().getProperty('max_connections') == '100'
        fs.properties().getProperty('max_error_retry') == '3'
        fs.properties().getProperty('upload_max_attempts') == '4'
        fs.properties().getProperty('upload_retry_sleep') == '200'
        fs.properties().getProperty('upload_chunk_size') == '7340032' //7MB
        fs.properties().getProperty('upload_max_threads') == '15'
        fs.properties().getProperty('max_concurrency') == '10'
        fs.properties().getProperty('max_native_memory') == '524288000' //500MB
        fs.properties().getProperty('minimum_part_size') == '7340032' //7MB
        fs.properties().getProperty('multipart_threshold') == '33554432' //32MB
    }

    def 'should resume an in-progress multipart upload'() {
        given:
        def awsClient = Mock(software.amazon.awssdk.services.s3.S3Client)
        def nfClient = Mock(nextflow.cloud.aws.nio.S3Client) {
            getClient() >> awsClient
        }
        def fs = Mock(S3FileSystem) {
            getClient() >> nfClient
            properties() >> new Properties()
        }
        def s3Path = new S3Path(fs, '/bucket/key.txt')
        def provider = new S3FileSystemProvider()

        def upload = MultipartUpload.builder().key('key.txt').uploadId('upload-id').initiated(Instant.now()).build()
        def part1 = Part.builder().partNumber(1).size(10).eTag('etag1').build()
        def part2 = Part.builder().partNumber(2).size(20).eTag('etag2').build()

        when:
        def handle = provider.resumeUpload(s3Path)

        then:
        1 * awsClient.listMultipartUploads(_) >> ListMultipartUploadsResponse.builder().uploads([upload]).build()
        1 * awsClient.listParts(_) >> ListPartsResponse.builder().parts([part1, part2]).build()

        and:
        handle != null
        handle.committedBytes() == 30
        handle.outputStream() instanceof S3OutputStream
    }

    def 'should paginate listParts when resuming a multipart upload'() {
        given:
        def awsClient = Mock(software.amazon.awssdk.services.s3.S3Client)
        def nfClient = Mock(nextflow.cloud.aws.nio.S3Client) {
            getClient() >> awsClient
        }
        def fs = Mock(S3FileSystem) {
            getClient() >> nfClient
            properties() >> new Properties()
        }
        def s3Path = new S3Path(fs, '/bucket/key.txt')
        def provider = new S3FileSystemProvider()

        def upload = MultipartUpload.builder().key('key.txt').uploadId('upload-id').initiated(Instant.now()).build()
        def part1 = Part.builder().partNumber(1).size(10).eTag('etag1').build()
        def part2 = Part.builder().partNumber(2).size(20).eTag('etag2').build()

        when:
        def handle = provider.resumeUpload(s3Path)

        then:
        1 * awsClient.listMultipartUploads(_) >> ListMultipartUploadsResponse.builder().uploads([upload]).build()
        2 * awsClient.listParts(_) >> { ListPartsRequest req ->
            req.partNumberMarker() == null
                ? ListPartsResponse.builder().parts([part1]).isTruncated(true).nextPartNumberMarker(1).build()
                : ListPartsResponse.builder().parts([part2]).isTruncated(false).build()
        }

        and:
        handle != null
        handle.committedBytes() == 30
    }

    def 'should paginate listMultipartUploads when resuming'() {
        given:
        def awsClient = Mock(software.amazon.awssdk.services.s3.S3Client)
        def nfClient = Mock(nextflow.cloud.aws.nio.S3Client) {
            getClient() >> awsClient
        }
        def fs = Mock(S3FileSystem) {
            getClient() >> nfClient
            properties() >> new Properties()
        }
        def s3Path = new S3Path(fs, '/bucket/key.txt')
        def provider = new S3FileSystemProvider()

        def upload = MultipartUpload.builder().key('key.txt').uploadId('upload-id').initiated(Instant.now()).build()
        def part1 = Part.builder().partNumber(1).size(10).eTag('etag1').build()

        when:
        def handle = provider.resumeUpload(s3Path)

        then:
        2 * awsClient.listMultipartUploads(_) >> { ListMultipartUploadsRequest req ->
            req.keyMarker() == null
                ? ListMultipartUploadsResponse.builder().uploads([]).isTruncated(true).nextKeyMarker('key.txt').build()
                : ListMultipartUploadsResponse.builder().uploads([upload]).isTruncated(false).build()
        }
        1 * awsClient.listParts(_) >> ListPartsResponse.builder().parts([part1]).build()

        and:
        handle != null
        handle.committedBytes() == 10
    }

    def 'should propagate S3 errors when resuming'() {
        given:
        def awsClient = Mock(software.amazon.awssdk.services.s3.S3Client)
        def nfClient = Mock(nextflow.cloud.aws.nio.S3Client) {
            getClient() >> awsClient
        }
        def fs = Mock(S3FileSystem) {
            getClient() >> nfClient
            properties() >> new Properties()
        }
        def s3Path = new S3Path(fs, '/bucket/key.txt')
        def provider = new S3FileSystemProvider()

        when:
        provider.resumeUpload(s3Path)

        then:
        1 * awsClient.listMultipartUploads(_) >> { throw S3Exception.builder().statusCode(503).build() }
        thrown(IOException)
    }
}
