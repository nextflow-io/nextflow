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

package nextflow.cloud.aws.nio;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.ByteBuffer;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.Phaser;
import java.util.concurrent.atomic.AtomicInteger;

import com.amazonaws.AmazonClientException;
import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.model.AbortMultipartUploadRequest;
import com.amazonaws.services.s3.model.CannedAccessControlList;
import com.amazonaws.services.s3.model.CompleteMultipartUploadRequest;
import com.amazonaws.services.s3.model.InitiateMultipartUploadRequest;
import com.amazonaws.services.s3.model.InitiateMultipartUploadResult;
import com.amazonaws.services.s3.model.ObjectMetadata;
import com.amazonaws.services.s3.model.ObjectTagging;
import com.amazonaws.services.s3.model.PartETag;
import com.amazonaws.services.s3.model.PutObjectRequest;
import com.amazonaws.services.s3.model.S3ObjectId;
import com.amazonaws.services.s3.model.SSEAlgorithm;
import com.amazonaws.services.s3.model.SSEAwsKeyManagementParams;
import com.amazonaws.services.s3.model.StorageClass;
import com.amazonaws.services.s3.model.Tag;
import com.amazonaws.services.s3.model.UploadPartRequest;
import com.amazonaws.util.Base64;
import nextflow.cloud.aws.nio.util.ByteBufferInputStream;
import nextflow.cloud.aws.nio.util.S3MultipartOptions;
import nextflow.util.Duration;
import nextflow.util.ThreadPoolHelper;
import nextflow.util.ThreadPoolManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import static java.util.Objects.requireNonNull;

/**
 * Parallel S3 multipart uploader. Based on the following code request
 * See https://github.com/Upplication/Amazon-S3-FileSystem-NIO2/pulls
 *
 * @author Paolo Di Tommaso
 * @author Tom Wieczorek
 */

public final class S3OutputStream extends OutputStream {


    private static final Logger log = LoggerFactory.getLogger(S3OutputStream.class);

    /**
     * Minimum multipart chunk size 5MB
     * https://docs.aws.amazon.com/AmazonS3/latest/userguide/qfacts.html
     */
    private static final int MIN_MULTIPART_UPLOAD = 5 * 1024 * 1024;

    /**
     * Amazon S3 API implementation to use.
     */
    private final AmazonS3 s3;

    /**
     * ID of the S3 object to store data into.
     */
    private final S3ObjectId objectId;

    /**
     * Amazon S3 storage class to apply to the newly created S3 object, if any.
     */
    private StorageClass storageClass;

    private SSEAlgorithm storageEncryption;

    private String kmsKeyId;

    private String contentType;

    /**
     * Indicates if the stream has been closed.
     */
    private volatile boolean closed;

    /**
     * Indicates if the upload has been aborted
     */
    private volatile boolean aborted;

    /**
     * If a multipart upload is in progress, holds the ID for it, {@code null} otherwise.
     */
    private volatile String uploadId;

    /**
     * If a multipart upload is in progress, holds the ETags of the uploaded parts, {@code null} otherwise.
     */
    private Queue<PartETag> partETags;

    /**
     * Holds upload request metadata
     */
    private final S3MultipartOptions request;

    /**
     * Instead of allocate a new buffer for each chunks recycle them, putting
     * a buffer instance into this queue when the upload process is completed
     */
    final private Queue<ByteBuffer> bufferPool = new ConcurrentLinkedQueue<ByteBuffer>();

    /**
     * The executor service (thread pool) which manages the upload in background
     */
    private ExecutorService executor;

    /**
     * The current working buffer
     */
    private ByteBuffer buf;

    private MessageDigest md5;

    /**
     * Phaser object to synchronize stream termination
     */
    private Phaser phaser;

    /**
     * Count the number of uploaded chunks
     */
    private int partsCount;

    private int bufferSize;

    private CannedAccessControlList cannedAcl;

    private List<Tag> tags;

    private final AtomicInteger bufferCounter = new AtomicInteger();

    /**
     * Creates a new {@code S3OutputStream} that writes data directly into the S3 object with the given {@code objectId}.
     * No special object metadata or storage class will be attached to the object.
     *
     */
    public S3OutputStream(final AmazonS3 s3, S3ObjectId objectId, S3MultipartOptions request) {
        this.s3 = requireNonNull(s3);
        this.objectId = requireNonNull(objectId);
        this.request = request;
        this.bufferSize = request.getBufferSize();
    }

    private ByteBuffer expandBuffer(ByteBuffer byteBuffer) {
        
        final float expandFactor = 2.5f;
        final int newCapacity = Math.min( (int)(byteBuffer.capacity() * expandFactor), bufferSize );

        // cast to prevent Java 8 / Java 11 cross compile-runtime error
        // https://www.morling.dev/blog/bytebuffer-and-the-dreaded-nosuchmethoderror/
        ((java.nio.Buffer)byteBuffer).flip();
        ByteBuffer expanded = ByteBuffer.allocate(newCapacity);
        expanded.order(byteBuffer.order());
        expanded.put(byteBuffer);
        return expanded;
    }

    public S3OutputStream setCannedAcl(CannedAccessControlList acl) {
        this.cannedAcl = acl;
        return this;
    }

    public S3OutputStream setTags(List<Tag> tags) {
        this.tags = tags;
        return this;
    }

    public S3OutputStream setStorageClass(String storageClass) {
        if( storageClass!=null )
            this.storageClass = StorageClass.fromValue(storageClass);
        return this;
    }

    public S3OutputStream setStorageEncryption(String storageEncryption) {
        if( storageEncryption!=null )
            this.storageEncryption = SSEAlgorithm.fromString(storageEncryption);
        return this;
    }

    public S3OutputStream setKmsKeyId(String kmsKeyId) {
        this.kmsKeyId = kmsKeyId;
        return this;
    }

    public S3OutputStream setContentType(String type) {
        this.contentType = type;
        return this;
    }

    /**
     * @return A MD5 message digester
     */
    private MessageDigest createMd5() {
        try {
            return MessageDigest.getInstance("MD5");
        }
        catch(NoSuchAlgorithmException e) {
            throw new IllegalStateException("Cannot find a MD5 algorithm provider",e);
        }
    }

    /**
     * Writes a byte into the uploader buffer. When it is full starts the upload process
     * in a asynchronous manner
     *
     * @param b The byte to be written
     * @throws IOException
     */
    @Override
    public void write (int b) throws IOException {
        if( closed ){
            throw new IOException("Can't write into a closed stream");
        }
        if( buf == null ) {
            buf = allocate();
            md5 = createMd5();
        }
        else if( !buf.hasRemaining() ) {
            if( buf.position() < bufferSize ) {
                buf = expandBuffer(buf);
            }
            else {
                flush();
                // create a new buffer
                buf = allocate();
                md5 = createMd5();
            }
        }

        buf.put((byte) b);
        // update the md5 checksum
        md5.update((byte) b);
    }

    /**
     * Flush the current buffer uploading to S3 storage
     *
     * @throws IOException
     */
    @Override
    public void flush() throws IOException {
        // send out the current buffer
        if( uploadBuffer(buf, false) ) {
            // clear the current buffer
            buf = null;
            md5 = null;
        }
    }

    private ByteBuffer allocate() {

        if( partsCount==0 ) {
            // this class is expected to be used to upload small files
            // start with a small buffer and growth if more space if necessary
            final int initialSize = 100 * 1024;
            return ByteBuffer.allocate(initialSize);
        }

        // try to reuse a buffer from the poll
        ByteBuffer result = bufferPool.poll();
        if( result != null ) {
            result.clear();
        }
        else {
            // allocate a new buffer
            log.debug("Allocating new buffer of {} bytes, total buffers {}", bufferSize, bufferCounter.incrementAndGet());
            result = ByteBuffer.allocate(bufferSize);
        }

        return result;
    }


    /**
     * Upload the given buffer to S3 storage in a asynchronous manner.
     * NOTE: when the executor service is busy (i.e. there are any more free threads)
     * this method will block
     *
     * return: true if the buffer can be reused, false if still needs to be used
     */
    private boolean uploadBuffer(ByteBuffer buf, boolean last) throws IOException {
        // when the buffer is empty nothing to do
        if( buf == null || buf.position()==0 ) { return false; }

        // Intermediate uploads needs to have at least MIN bytes
        if( buf.position() < MIN_MULTIPART_UPLOAD && !last){
            return false;
        }

        if (partsCount == 0) {
            init();
        }

        // set the buffer in read mode and submit for upload
        executor.submit( task(buf, md5.digest(), ++partsCount) );

        return true;
    }

    /**
     * Initialize multipart upload data structures
     *
     * @throws IOException
     */
    private void init() throws IOException {
        // get the upload id
        uploadId = initiateMultipartUpload().getUploadId();
        if (uploadId == null) {
            throw new IOException("Failed to get a valid multipart upload ID from Amazon S3");
        }
        // create the executor
        executor = getOrCreateExecutor(request.getMaxThreads());
        partETags = new LinkedBlockingQueue<>();
        phaser = new Phaser();
        phaser.register();
        log.trace("[S3 phaser] Register - Starting S3 upload: {}; chunk-size: {}; max-threads: {}", uploadId, bufferSize, request.getMaxThreads());
    }


    /**
     * Creates a {@link Runnable} task to handle the upload process
     * in background
     *
     * @param buffer The buffer to be uploaded
     * @param partIndex The index count
     * @return
     */
    private Runnable task(final ByteBuffer buffer, final byte[] checksum, final int partIndex) {

        phaser.register();
        log.trace("[S3 phaser] Task register");
        return new Runnable() {
            @Override
            public void run() {
                try {
                    uploadPart(buffer, checksum, partIndex, false);
                }
                catch (IOException e) {
                    final StringWriter writer = new StringWriter();
                    e.printStackTrace(new PrintWriter(writer));
                    log.error("Upload: {} > Error for part: {}\nCaused by: {}", uploadId, partIndex, writer.toString());
                }
                finally {
                    log.trace("[S3 phaser] Task arriveAndDeregisterphaser");
                    phaser.arriveAndDeregister();
                }
            }
        };

    }

    /**
     * Close the stream uploading any remaining buffered data
     *
     * @throws IOException
     */
    @Override
    public void close() throws IOException {
        if (closed) {
            return;
        }

        if (uploadId == null) {
            if( buf != null )
                putObject(buf, md5.digest());
            else
                // this is needed when trying to upload an empty 
                putObject(new ByteArrayInputStream(new byte[]{}), 0, createMd5().digest());
        }
        else {
            // -- upload remaining chunk
            if( buf != null )
                uploadBuffer(buf, true);

            // -- shutdown upload executor and await termination
            log.trace("[S3 phaser] Close arriveAndAwaitAdvance");
            phaser.arriveAndAwaitAdvance();

            // -- complete upload process
            completeMultipartUpload();
        }

        closed = true;
    }

    /**
     * Starts the multipart upload process
     *
     * @return An instance of {@link InitiateMultipartUploadResult}
     * @throws IOException
     */
    private InitiateMultipartUploadResult initiateMultipartUpload() throws IOException {
        final InitiateMultipartUploadRequest request = //
                new InitiateMultipartUploadRequest(objectId.getBucket(), objectId.getKey());
        final ObjectMetadata metadata = new ObjectMetadata();

        if (storageClass != null) {
            request.setStorageClass(storageClass);
        }

        if( cannedAcl != null ) {
            request.withCannedACL(cannedAcl);
        }

        if( kmsKeyId !=null ) {
            request.withSSEAwsKeyManagementParams( new SSEAwsKeyManagementParams(kmsKeyId) );
        }

        if( storageEncryption != null ) {
            metadata.setSSEAlgorithm(storageEncryption.toString());
            request.setObjectMetadata(metadata);
        }

        if( contentType != null ) {
            metadata.setContentType(contentType);
            request.setObjectMetadata(metadata);
        }

        if( log.isTraceEnabled() ) {
            log.trace("S3 initiateMultipartUpload {}", request);
        }

        try {
            return s3.initiateMultipartUpload(request);
        } catch (final AmazonClientException e) {
            throw new IOException("Failed to initiate Amazon S3 multipart upload", e);
        }
    }

    /**
     * Upload the given buffer to the S3 storage using a multipart process
     *
     * @param buf The buffer holding the data to upload
     * @param partNumber The progressive index of this chunk (1-based)
     * @param lastPart {@code true} when it is the last chunk
     * @throws IOException
     */
    private void uploadPart( final ByteBuffer buf, final byte[] checksum, final int partNumber, final boolean lastPart ) throws IOException {
        // cast to prevent Java 8 / Java 11 cross compile-runtime error
        // https://www.morling.dev/blog/bytebuffer-and-the-dreaded-nosuchmethoderror/
        ((java.nio.Buffer)buf).flip();
        ((java.nio.Buffer)buf).mark();

        int attempt=0;
        boolean success=false;
        try {
            while( !success ) {
                attempt++;
                int len = buf.limit();
                try {
                    log.trace("Uploading part {} with length {} attempt {} for {} ", partNumber, len, attempt, objectId);
                    uploadPart( new ByteBufferInputStream(buf), len, checksum , partNumber, lastPart );
                    success=true;
                }
                catch (AmazonClientException | IOException e) {
                    if( attempt == request.getMaxAttempts() )
                        throw new IOException("Failed to upload multipart data to Amazon S3", e);

                    log.debug("Failed to upload part {} attempt {} for {} -- Caused by: {}", partNumber, attempt, objectId, e.getMessage());
                    sleep(request.getRetrySleep());
                    buf.reset();
                }
            }
        }
        finally {
            if (!success) {
                closed = true;
                abortMultipartUpload();
            }
            bufferPool.offer(buf);
        }

    }

    private void uploadPart(final InputStream content, final long contentLength, final byte[] checksum, final int partNumber, final boolean lastPart)
            throws IOException {

        if (aborted) return;

        final UploadPartRequest request = new UploadPartRequest();
        request.setBucketName(objectId.getBucket());
        request.setKey(objectId.getKey());
        request.setUploadId(uploadId);
        request.setPartNumber(partNumber);
        request.setPartSize(contentLength);
        request.setInputStream(content);
        request.setLastPart(lastPart);
        request.setMd5Digest(Base64.encodeAsString(checksum));

        final PartETag partETag = s3.uploadPart(request).getPartETag();
        log.trace("Uploaded part {} with length {} for {}: {}", partETag.getPartNumber(), contentLength, objectId, partETag.getETag());
        partETags.add(partETag);

    }

    private void sleep( long millis ) {
        try {
            Thread.sleep(millis);
        }
        catch (InterruptedException e) {
            log.trace("Sleep was interrupted -- Cause: {}", e.getMessage());
        }
    }

    /**
     * Aborts the multipart upload process
     */
    private synchronized void abortMultipartUpload() {
        if (aborted) return;

        log.debug("Aborting multipart upload {} for {}", uploadId, objectId);
        try {
            s3.abortMultipartUpload(new AbortMultipartUploadRequest(objectId.getBucket(), objectId.getKey(), uploadId));
        }
        catch (final AmazonClientException e) {
            log.warn("Failed to abort multipart upload {}: {}", uploadId, e.getMessage());
        }
        aborted = true;
        log.trace("[S3 phaser] MultipartUpload arriveAndDeregister");
        phaser.arriveAndDeregister();
    }

    /**
     * Completes the multipart upload process
     * @throws IOException
     */
    private void completeMultipartUpload() throws IOException {
        // if aborted upload just ignore it
        if( aborted ) return;

        final int partCount = partETags.size();
        log.trace("Completing upload to {} consisting of {} parts", objectId, partCount);

        try {
            s3.completeMultipartUpload(new CompleteMultipartUploadRequest( //
                    objectId.getBucket(), objectId.getKey(), uploadId, new ArrayList<>(partETags)));
        } catch (final AmazonClientException e) {
            throw new IOException("Failed to complete Amazon S3 multipart upload", e);
        }

        log.trace("Completed upload to {} consisting of {} parts", objectId, partCount);

        uploadId = null;
        partETags = null;
    }

    /**
     * Stores the given buffer using a single-part upload process
     * @param buf
     * @throws IOException
     */
    private void putObject(ByteBuffer buf, byte[] checksum) throws IOException {
        // cast to prevent Java 8 / Java 11 cross compile-runtime error
        // https://www.morling.dev/blog/bytebuffer-and-the-dreaded-nosuchmethoderror/
        ((java.nio.Buffer)buf).flip();
        putObject(new ByteBufferInputStream(buf), buf.limit(), checksum);
    }

    /**
     * Stores the given buffer using a single-part upload process
     *
     * @param contentLength
     * @param content
     * @throws IOException
     */
    private void putObject(final InputStream content, final long contentLength, byte[] checksum) throws IOException {

        final ObjectMetadata meta = new ObjectMetadata();
        meta.setContentLength(contentLength);
        meta.setContentMD5( Base64.encodeAsString(checksum) );

        final PutObjectRequest request = new PutObjectRequest(objectId.getBucket(), objectId.getKey(), content, meta);
        if( cannedAcl!=null ) {
            request.withCannedAcl(cannedAcl);
        }

        if (storageClass != null) {
            request.setStorageClass(storageClass);
        }

        if( tags!=null && tags.size()>0 ) {
            request.setTagging( new ObjectTagging(tags) );
        }

        if( kmsKeyId !=null ) {
            request.withSSEAwsKeyManagementParams( new SSEAwsKeyManagementParams(kmsKeyId) );
        }

        if( storageEncryption != null ) {
            meta.setSSEAlgorithm( storageEncryption.toString() );
        }

        if( contentType != null ) {
            meta.setContentType(contentType);
        }

        if( log.isTraceEnabled() ) {
            log.trace("S3 putObject {}", request);
        }

        try {
            s3.putObject(request);
        } catch (final AmazonClientException e) {
            throw new IOException("Failed to put data into Amazon S3 object", e);
        }
    }

    /**
     * @return Number of uploaded chunks
     */
    int getPartsCount() {
        return partsCount;
    }


    /** holds a singleton executor instance */
    static private volatile ExecutorService executorSingleton;

    /**
     * Creates a singleton executor instance.
     *
     * @param maxThreads
     *          The max number of allowed threads in the executor pool.
     *          NOTE: changing the size parameter after the first invocation has no effect.
     * @return The executor instance
     */
    static synchronized ExecutorService getOrCreateExecutor(int maxThreads) {
        if( executorSingleton == null ) {
            executorSingleton = ThreadPoolManager.create("S3StreamUploader", maxThreads);
        }
        return executorSingleton;
    }

}
