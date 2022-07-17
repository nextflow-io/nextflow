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

package com.upplication.s3fs.util;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.List;

import com.amazonaws.AmazonClientException;
import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.model.CannedAccessControlList;
import com.amazonaws.services.s3.model.ObjectMetadata;
import com.amazonaws.services.s3.model.ObjectTagging;
import com.amazonaws.services.s3.model.PutObjectRequest;
import com.amazonaws.services.s3.model.S3ObjectId;
import com.amazonaws.services.s3.model.SSEAlgorithm;
import com.amazonaws.services.s3.model.SSEAwsKeyManagementParams;
import com.amazonaws.services.s3.model.StorageClass;
import com.amazonaws.services.s3.model.Tag;
import com.amazonaws.util.Base64;
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

public final class S3CopyStream extends OutputStream {

    private static final Logger log = LoggerFactory.getLogger(S3CopyStream.class);

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

    /**
     * Indicates if the stream has been closed.
     */
    private volatile boolean closed;

    /**
     * Indicates if the upload has been aborted
     */
    private volatile boolean aborted;

    private MessageDigest md5;

    private CannedAccessControlList cannedAcl;

    private List<Tag> tags;

    private CopyOutputStream buffer;

    /**
     * Creates a new {@code S3OutputStream} that writes data directly into the S3 object with the given {@code objectId}.
     * No special object metadata or storage class will be attached to the object.
     *
     */
    public S3CopyStream(final AmazonS3 s3, S3ObjectId objectId) {
        this.s3 = requireNonNull(s3);
        this.objectId = requireNonNull(objectId);
        this.md5 = createMd5();
        this.buffer = new CopyOutputStream();
    }

    public S3CopyStream setCannedAcl(CannedAccessControlList acl) {
        this.cannedAcl = acl;
        return this;
    }

    public S3CopyStream setTags(List<Tag> tags) {
        this.tags = tags;
        return this;
    }

    public S3CopyStream setStorageClass(String storageClass) {
        if( storageClass!=null )
            this.storageClass = StorageClass.fromValue(storageClass);
        return this;
    }

    public S3CopyStream setStorageEncryption(String storageEncryption) {
        if( storageEncryption!=null )
            this.storageEncryption = SSEAlgorithm.fromString(storageEncryption);
        return this;
    }

    public S3CopyStream setKmsKeyId(String kmsKeyId) {
        this.kmsKeyId = kmsKeyId;
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

    public void write(byte b[], int off, int len) throws IOException {
        if( closed ){
            throw new IOException("Can't write into a closed stream");
        }
        buffer.write(b,off,len);
        md5.update(b,off,len);
    }

    /**
     * Writes a byte into the uploader buffer. When it is full starts the upload process
     * in a asynchornous manner
     *
     * @param b The byte to be written
     * @throws IOException
     */
    @Override
    public void write (int b) throws IOException {
        if( closed ){
            throw new IOException("Can't write into a closed stream");
        }
        buffer.write((byte) b);
        md5.update((byte) b);
    }


    /**
     * Close the stream uploading any remaning buffered data
     *
     * @throws IOException
     */
    @Override
    public void close() throws IOException {
        if (closed) {
            return;
        }

        putObject(buffer.toInputStream(), buffer.size(), md5.digest());
        closed = true;
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

        if( log.isTraceEnabled() ) {
            log.trace("S3 putObject {}", request);
        }

        try {
            s3.putObject(request);
        } catch (final AmazonClientException e) {
            throw new IOException("Failed to put data into Amazon S3 object", e);
        }
    }

}
