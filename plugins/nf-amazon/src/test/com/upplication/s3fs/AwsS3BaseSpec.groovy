package com.upplication.s3fs

import java.nio.ByteBuffer
import java.nio.channels.SeekableByteChannel
import java.nio.file.Path
import java.nio.file.Paths

import com.amazonaws.services.s3.AmazonS3
import com.amazonaws.services.s3.model.AmazonS3Exception
import com.amazonaws.services.s3.model.ListVersionsRequest
import com.amazonaws.services.s3.model.S3ObjectSummary
import com.amazonaws.services.s3.model.S3VersionSummary
import org.slf4j.Logger
import org.slf4j.LoggerFactory
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
trait AwsS3BaseSpec {

    static final Logger log = LoggerFactory.getLogger(AwsS3BaseSpec)

    abstract AmazonS3 getS3Client()

    S3Path s3path(String path) {
        path = path.replaceAll(~'^s3://(?!/)','s3:///')
        return (S3Path) Path.of(new URI(path))
    }

    String createBucket(String bucketName) {
        s3Client.createBucket(bucketName)
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
        return s3Client .putObject(bucketName, blobName, content)
    }

    def createDirectory(String path) {
        log.debug "Creating blob directory '$path'"
        def (bucketName, blobName) = splitName(path)
        blobName += '/'
        s3Client.putObject(bucketName, blobName, '')
    }

    def deleteObject(String path) {
        log.debug "Deleting blob object '$path'"
        def (bucketName, blobName) = splitName(path)
        blobName += '/'
        s3Client.deleteObject(bucketName, blobName)
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
        def objectListing = s3Client.listObjects(bucketName);
        while (true) {
            Iterator<S3ObjectSummary> objIter = objectListing.getObjectSummaries().iterator();
            while (objIter.hasNext()) {
                s3Client.deleteObject(bucketName, objIter.next().getKey());
            }

            // If the bucket contains many objects, the listObjects() call
            // might not return all of the objects in the first listing. Check to
            // see whether the listing was truncated. If so, retrieve the next page of objects
            // and delete them.
            if (objectListing.isTruncated()) {
                objectListing = s3Client.listNextBatchOfObjects(objectListing);
            } else {
                break;
            }
        }

        // Delete all object versions (required for versioned buckets).
        def versionList = s3Client.listVersions(new ListVersionsRequest().withBucketName(bucketName));
        while (true) {
            Iterator<S3VersionSummary> versionIter = versionList.getVersionSummaries().iterator();
            while (versionIter.hasNext()) {
                S3VersionSummary vs = versionIter.next();
                s3Client.deleteVersion(bucketName, vs.getKey(), vs.getVersionId());
            }

            if (versionList.isTruncated()) {
                versionList = s3Client.listNextBatchOfVersions(versionList);
            } else {
                break;
            }
        }

        // After all objects and object versions are deleted, delete the bucket.
        s3Client.deleteBucket(bucketName);

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
                return s3Client.doesBucketExist(path.getName(0).toString())
            }
            else {
                s3Client.getObject(bucketName, blobName).getObjectMetadata()
                return true
            }
        }
        catch (AmazonS3Exception e) {
            if( e.statusCode == 404 )
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
                .getObject(bucketName, blobName)
                .getObjectContent()
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
