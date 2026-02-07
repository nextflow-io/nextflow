/*
 * Copyright 2020-2025, Seqera Labs
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

package nextflow.cloud.aws.nio

import nextflow.Global
import nextflow.Session
import nextflow.cloud.aws.nio.util.S3MultipartOptions
import nextflow.file.FileHelper
import software.amazon.awssdk.services.s3.S3Client
import software.amazon.awssdk.services.s3.model.CompleteMultipartUploadRequest
import software.amazon.awssdk.services.s3.model.CreateMultipartUploadResponse
import software.amazon.awssdk.services.s3.model.UploadPartResponse
import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Specification

import java.nio.file.Files
import java.nio.file.attribute.BasicFileAttributes

/**
 * Test for S3OutputStream
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class S3OutputStreamTest extends Specification implements AwsS3BaseSpec {

    private S3Client s3Client0

    S3Client getS3Client() { s3Client0 }

    static private Map config0() {
        def accessKey = System.getenv('AWS_S3FS_ACCESS_KEY')
        def secretKey = System.getenv('AWS_S3FS_SECRET_KEY')
        return [aws: [accessKey: accessKey, secretKey: secretKey]]
    }

    def setup() {
        def fs = (S3FileSystem) FileHelper.getOrCreateFileSystemFor(URI.create("s3:///"), config0().aws)
        s3Client0 = fs.client.getClient()
        and:
        def cfg = config0()
        Global.config = cfg
        Global.session = Mock(Session) { getConfig() >> cfg }
    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    @Requires({System.getenv('AWS_S3FS_ACCESS_KEY') && System.getenv('AWS_S3FS_SECRET_KEY')})
    def 'should ensure multipart is used'() {
        given:
        def bucket = createBucket()
        and:
        def chunksize = 6 * 1024 * 1024
        def bytes = new byte[chunksize]
        new Random().nextBytes(bytes)
        final path = s3path("s3://$bucket/file.txt")
        def multipart = new S3MultipartOptions()
        multipart.setChunkSize(chunksize)
        multipart.setBufferSize(chunksize)
        when:
        def writer = new S3OutputStream(s3Client0, path.toS3ObjectId(), multipart)
        10.times { it ->
            writer.write(bytes);
            writer.flush()
        }
        writer.close()

        then:
        writer.partsCount == 10
        existsPath(path)
        Files.readAttributes(path, BasicFileAttributes).size() == 10 * chunksize

        cleanup:
        if( bucket ) deleteBucket(bucket)
    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    @Requires({System.getenv('AWS_S3FS_ACCESS_KEY') && System.getenv('AWS_S3FS_SECRET_KEY')})
    def 'should upload empty stream'() {
        given:
        def bucket = createBucket()
        and:
        final path = s3path("s3://$bucket/file.txt")
        def multipart = new S3MultipartOptions()
        when:
        def writer = new S3OutputStream(s3Client0, path.toS3ObjectId(), multipart)
        writer.close()

        then:
        writer.partsCount == 0
        existsPath(path)
        Files.readAttributes(path, BasicFileAttributes).size() == 0

        cleanup:
        if( bucket ) deleteBucket(bucket)
    }
    @IgnoreIf({System.getenv('NXF_SMOKE')})
    @Requires({System.getenv('AWS_S3FS_ACCESS_KEY') && System.getenv('AWS_S3FS_SECRET_KEY')})
    def 'should upload without multipart'() {
        given:
        def bucket = createBucket()
        and:
        def TEXT = randomText(50 * 1024)
        final path = s3path("s3://$bucket/file.txt")
        def multipart = new S3MultipartOptions()
        when:
        def writer = new S3OutputStream(s3Client0, path.toS3ObjectId(), multipart)
        writer.write(TEXT.bytes)
        writer.close()

        then:
        writer.partsCount == 0
        existsPath(path)
        readObject(path) == TEXT

        cleanup:
        if( bucket ) deleteBucket(bucket)
    }

    def 'should send sorted parts to completeMultipartUpload'() {
        given:
        final path = s3path("s3://test/file.txt")
        def multipart = new S3MultipartOptions()
        def client = Mock(S3Client)
        def capturedParts = null

        def writer = new S3OutputStream(client, path.toS3ObjectId(), multipart)

        when: 'simulate unsorted uploads'
        writer.init()
        writer.uploadPart(InputStream.nullInputStream(), 25, "checksum".bytes, 2, true)
        writer.uploadPart(InputStream.nullInputStream(), 25, "checksum".bytes, 0, false)
        writer.uploadPart(InputStream.nullInputStream(), 25, "checksum".bytes, 1, false)
        writer.completeMultipartUpload()

        then:
        1 * client.createMultipartUpload(_) >> CreateMultipartUploadResponse.builder().uploadId("upload-id").build()
        3 * client.uploadPart(_,_) >> { UploadPartResponse.builder().eTag('etag').build()}
        1 * client.completeMultipartUpload(_ as CompleteMultipartUploadRequest) >> { CompleteMultipartUploadRequest req ->
            capturedParts = req.multipartUpload().parts()
            return null
        }
        capturedParts[0].partNumber() == 0
        capturedParts[1].partNumber() == 1
        capturedParts[2].partNumber() == 2

    }

    def 'should include CRC32C checksum for uploads'() {
        given:
        def path = s3path("s3://test-bucket--use1-az1--x-s3/file.txt")
        def multipart = new S3MultipartOptions()
        def client = Mock(S3Client)
        def writer = new S3OutputStream(client, path.toS3ObjectId(), multipart)
        def capturedUploadRequest = null
        def capturedPutRequest = null
        // Pre-computed CRC32C checksum bytes (4 bytes big-endian)
        def checksumBytes = [0x12, 0x34, 0x56, 0x78] as byte[]
        def expectedChecksum = Base64.getEncoder().encodeToString(checksumBytes)

        when: 'upload with multipart'
        writer.init()
        writer.uploadPart(InputStream.nullInputStream(), 25, checksumBytes, 1, true)

        then:
        1 * client.createMultipartUpload(_) >> CreateMultipartUploadResponse.builder().uploadId("upload-id").build()
        1 * client.uploadPart(_, _) >> { args ->
            capturedUploadRequest = args[0]
            return UploadPartResponse.builder().eTag('etag').build()
        }
        capturedUploadRequest.contentMD5() == null
        capturedUploadRequest.checksumCRC32C() == expectedChecksum

        when: 'upload without multipart'
        path = s3path("s3://test-bucket/small.txt")
        writer = new S3OutputStream(client, path.toS3ObjectId(), multipart)
        writer.write("test".bytes)
        writer.close()

        then:
        1 * client.putObject(_, _) >> { args ->
            capturedPutRequest = args[0]
            return null
        }
        capturedPutRequest.contentMD5() == null
        capturedPutRequest.checksumCRC32C() != null
        and: 'verify CRC32C is correct for "test" bytes'
        def crc = new java.util.zip.CRC32C()
        crc.update("test".bytes)
        def expectedCrc = crc.getValue()
        def expectedCrcBytes = [
            (byte)((expectedCrc >> 24) & 0xFF),
            (byte)((expectedCrc >> 16) & 0xFF),
            (byte)((expectedCrc >> 8) & 0xFF),
            (byte)(expectedCrc & 0xFF)
        ] as byte[]
        capturedPutRequest.checksumCRC32C() == Base64.getEncoder().encodeToString(expectedCrcBytes)
    }
}
