package nextflow.cloud.azure.nio

import java.nio.ByteBuffer
import java.nio.channels.SeekableByteChannel
import java.nio.file.Path
import java.nio.file.Paths

import com.azure.storage.blob.BlobServiceClient
import com.azure.storage.blob.models.BlobStorageException
import org.slf4j.Logger
import org.slf4j.LoggerFactory
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
trait AzBaseSpec {

    static final Logger log = LoggerFactory.getLogger(AzBaseSpec)

    abstract BlobServiceClient getStorageClient()

    String createBucket(String containerName) {
        storageClient.createBlobContainer(containerName)
        return containerName
    }

    String createBucket() {
        def name = getRndBucketName()
        log.debug "Creating blob bucket '$name'"
        createBucket(name)
    }

    String getRndBucketName() {
        return "nf-az-blob-${UUID.randomUUID()}"
    }

    def createObject(String path, String content) {
        createObject(Paths.get(path), content)
    }

    def createObject(Path path, String content) {
        log.debug "Creating blob object '$path'"
        if( path.nameCount<=1 )
            throw new IllegalArgumentException("There should be at least one dir level: $path")
        final containerName = path.subpath(0,1).toString()
        final blobName = path.subpath(1, path.nameCount).toString()
        final blob = storageClient
                .getBlobContainerClient(containerName)
                .getBlobClient(blobName)
                .getAppendBlobClient()

        blob.create()

        blob
                .blobOutputStream
                .withStream { it.write(content.bytes) }
    }

    def createDirectory(String path) {
        log.debug "Creating blob directory '$path'"
        final containerName = path.tokenize('/')[0]
        final blobName = path.tokenize('/')[1..-1].join('/') + '/'
        storageClient
                .getBlobContainerClient(containerName)
                .getBlobClient(blobName)
                .appendBlobClient
                .create()
    }

    def deleteObject(String path) {
        log.debug "Deleting blob object '$path'"
        final containerName = path.tokenize('/')[0]
        final blobName = path.tokenize('/')[1..-1].join('/') + '/'

        storageClient
                .getBlobContainerClient(containerName)
                .getBlobClient(blobName)
                .delete()
    }

    def deleteBucket(Path path) {
        log.debug "Deleting blob bucket '$path'"
        assert path.nameCount == 1
        deleteBucket(path.getName(0).toString())
    }

    def deleteBucket(String bucketName) {
        log.debug "Deleting blob bucket '$bucketName'"
        if( !bucketName )
            return

        storageClient.deleteBlobContainer(bucketName)
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
        if( path.nameCount == 1 ) {
            try {
                storageClient
                        .getBlobContainerClient(path.getName(0).toString())
                        .getProperties()
                return true
            }
            catch (BlobStorageException e) {
                if( e.statusCode == 404 )
                    return false
                throw e
            }
        }

        final containerName = path.subpath(0,1).toString()
        final blobName = path.subpath(1,path.nameCount).toString()
        try {
            storageClient.getBlobContainerClient(containerName)
                .getBlobClient(blobName)
                .getProperties()
            return true
        }
        catch (BlobStorageException e) {
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
        final containerName = path.subpath(0,1).toString()
        final blobName = path.subpath(1,path.nameCount).toString()
        try {
            return storageClient.getBlobContainerClient(containerName)
                    .getBlobClient(blobName)
                    .openInputStream()
                    .text

        }
        catch (BlobStorageException e) {
            throw e
        }
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
