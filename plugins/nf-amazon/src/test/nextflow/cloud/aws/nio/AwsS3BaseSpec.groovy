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

import software.amazon.awssdk.core.sync.RequestBody
import software.amazon.awssdk.services.s3.model.GetObjectRequest
import software.amazon.awssdk.services.s3.model.HeadBucketRequest
import software.amazon.awssdk.services.s3.model.HeadObjectRequest

import java.nio.ByteBuffer
import java.nio.channels.SeekableByteChannel
import java.nio.file.Path
import java.nio.file.Paths

import software.amazon.awssdk.services.s3.S3Client
import software.amazon.awssdk.services.s3.model.CreateBucketRequest
import software.amazon.awssdk.services.s3.model.DeleteBucketRequest
import software.amazon.awssdk.services.s3.model.DeleteObjectRequest
import software.amazon.awssdk.services.s3.model.S3Exception
import software.amazon.awssdk.services.s3.model.ListObjectsV2Request
import software.amazon.awssdk.services.s3.model.ListObjectVersionsRequest
import software.amazon.awssdk.services.s3.model.S3Object
import software.amazon.awssdk.services.s3.model.ObjectVersion
import software.amazon.awssdk.services.s3.model.PutObjectRequest
import nextflow.cloud.aws.util.S3PathFactory
import org.slf4j.Logger
import org.slf4j.LoggerFactory
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
trait AwsS3BaseSpec {

    static final Logger log = LoggerFactory.getLogger(AwsS3BaseSpec)

    abstract S3Client getS3Client()

    S3Path s3path(String path) {
        return (S3Path) S3PathFactory.parse(path)
    }

    String createBucket(String bucketName) {
        s3Client.createBucket(CreateBucketRequest.builder().bucket(bucketName).build() as CreateBucketRequest)
        return bucketName
    }

    String createBucket() {
        def name = getRndBucketName()
        log.debug "Creating s3 bucket '$name'"
        createBucket(name)
    }

    String getRndBucketName() {
        return "nf-s3fs-test-${UUID.randomUUID()}"
    }

    def createObject(String path, String content) {
        createObject(Paths.get(path), content)
    }

    private List<String> splitName(path) {
        def items = path.toString().tokenize('/')
        return items.size()==1
            ? [ items[0], null ]
            : [ items[0], items[1..-1].join('/') ]
    }

    def createObject(Path path, String content) {
        log.debug "Creating s3 blob object '$path'"
        def (bucketName, blobName) = splitName(path)
        if( !blobName )
            throw new IllegalArgumentException("There should be at least one dir level: $path")
        return s3Client.putObject(PutObjectRequest.builder().bucket(bucketName).key(blobName).build() as PutObjectRequest, RequestBody.fromBytes(content.bytes))
    }

    def createDirectory(String path) {
        log.debug "Creating blob directory '$path'"
        def (bucketName, blobName) = splitName(path)
        blobName += '/'
        s3Client.putObject(PutObjectRequest.builder().bucket(bucketName).key(blobName).build() as PutObjectRequest, RequestBody.empty())
    }

    def deleteObject(String path) {
        log.debug "Deleting blob object '$path'"
        def (bucketName, blobName) = splitName(path)
        blobName += '/'
        s3Client.deleteObject(DeleteObjectRequest.builder().bucket(bucketName).key(blobName).build() as DeleteObjectRequest)
    }

    def deleteBucket(Path path) {
        log.debug "Deleting blob bucket '$path'"
        def (bucketName, blobName) = splitName(path)
        assert blobName == null
        deleteBucket(bucketName)
    }

    def deleteBucket(String bucketName) {
        log.debug "Deleting blob bucket '$bucketName'"
        if( !bucketName )
            return

        // Delete all objects from the bucket. This is sufficient
        // for unversioned buckets. For versioned buckets, when you attempt to delete objects, Amazon S3 inserts
        // delete markers for all objects, but doesn't delete the object versions.
        // To delete objects from versioned buckets, delete all of the object versions before deleting
        // the bucket (see below for an example).
        def objectListingIterator = s3Client.listObjectsV2Paginator(ListObjectsV2Request.builder().bucket(bucketName).build() as ListObjectsV2Request).iterator();
        while (objectListingIterator.hasNext()) {
            Iterator<S3Object> objIter = objectListingIterator.next().contents().iterator();
            while (objIter.hasNext()) {
                s3Client.deleteObject(DeleteObjectRequest.builder().bucket(bucketName).key(objIter.next().key()).build() as DeleteObjectRequest);
            }
        }

        // Delete all object versions (required for versioned buckets).
        def versionListIterator = s3Client.listObjectVersionsPaginator(ListObjectVersionsRequest.builder().bucket(bucketName).build() as ListObjectVersionsRequest).iterator();
        while ( versionListIterator.hasNext()){
            Iterator<ObjectVersion> versionIter = versionListIterator.next().versions().iterator();
            while ( versionIter.hasNext() ) {
                ObjectVersion vs = versionIter.next();
                s3Client.deleteObject(DeleteObjectRequest.builder().bucket(bucketName).key(vs.key()).versionId(vs.versionId()).build() as DeleteObjectRequest);
            }

        }

        // After all objects and object versions are deleted, delete the bucket.
        s3Client.deleteBucket( DeleteBucketRequest.builder().bucket(bucketName).build() as DeleteBucketRequest);

    }

    def tryDeleteBucket(String bucketName) {
        try {
            deleteBucket(bucketName)
        }
        catch (Throwable t) {
            log.warn ("Unable to delete blob bucket '$bucketName' - Raeason: ${t.message ?: t}")
        }
    }

    boolean existsPath(String path) {
        log.debug "Check blob path exists '$path'"
        existsPath(Paths.get(path))
    }

    boolean existsPath(Path path) {
        log.debug "Check blob path exists '$path'"
        def (bucketName, blobName) = splitName(path)
        if( !bucketName )
            throw new IllegalArgumentException("Invalid S3 path $path")

        try {
            if( !blobName ) {
                return s3Client.headBucket(HeadBucketRequest.builder().bucket(bucketName).build() as HeadBucketRequest)
            }
            else {
                s3Client.headObject(HeadObjectRequest.builder().bucket(bucketName).key(blobName).build() as HeadObjectRequest)
                return true
            }
        }
        catch (S3Exception e) {
            if( e.statusCode() == 404 )
                return false
            throw e
        }
    }

    String readObject(String path) {
        log.debug "Reading blob object '$path'"
        readObject(Paths.get(path))
    }

    String readObject(Path path) {
        log.debug "Reading blob object '$path'"
        def (bucketName, blobName) = splitName(path)
        return s3Client
                .getObject(GetObjectRequest.builder().bucket(bucketName).key(blobName).build() as GetObjectRequest)
                .getText()
    }


    String randomText(int size) {
        def result = new StringBuilder()
        while( result.size() < size ) {
            result << UUID.randomUUID().toString() << '\n'
        }
        return result.toString()
    }

    String readChannel(SeekableByteChannel sbc, int buffLen )  {
        def buffer = new ByteArrayOutputStream()
        ByteBuffer bf = ByteBuffer.allocate(buffLen)
        while((sbc.read(bf))>0) {
            bf.flip();
            buffer.write(bf.array(), 0, bf.limit())
            bf.clear();
        }

        buffer.toString()
    }

    void writeChannel( SeekableByteChannel channel, String content, int buffLen ) {

        def bytes = content.getBytes()
        ByteBuffer buf = ByteBuffer.allocate(buffLen);
        int i=0
        while( i < bytes.size()) {

            def len = Math.min(buffLen, bytes.size()-i);
            buf.clear();
            buf.put(bytes, i, len);
            buf.flip();
            channel.write(buf);

            i += len
        }

    }

}
