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
}
