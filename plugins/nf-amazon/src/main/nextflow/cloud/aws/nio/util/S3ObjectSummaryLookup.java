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

package nextflow.cloud.aws.nio.util;

import java.io.IOException;
import java.nio.file.NoSuchFileException;
import java.util.List;

import com.amazonaws.services.s3.model.AmazonS3Exception;
import com.amazonaws.services.s3.model.ListObjectsRequest;
import com.amazonaws.services.s3.model.ObjectListing;
import com.amazonaws.services.s3.model.ObjectMetadata;
import com.amazonaws.services.s3.model.S3Object;
import com.amazonaws.services.s3.model.S3ObjectSummary;
import nextflow.cloud.aws.nio.S3Client;
import nextflow.cloud.aws.nio.S3Path;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class S3ObjectSummaryLookup {

    private static final Logger log = LoggerFactory.getLogger(S3ObjectSummary.class);

    /**
     * Get the {@link com.amazonaws.services.s3.model.S3ObjectSummary} that represent this Path or her first child if this path not exists
     * @param s3Path {@link S3Path}
     * @return {@link com.amazonaws.services.s3.model.S3ObjectSummary}
     * @throws java.nio.file.NoSuchFileException if not found the path and any child
     */
    public S3ObjectSummary lookup(S3Path s3Path) throws NoSuchFileException {

        /*
         * check is object summary has been cached
         */
        S3ObjectSummary summary = s3Path.fetchObjectSummary();
        if( summary != null ) {
            return summary;
        }

        final S3Client client = s3Path.getFileSystem().getClient();

        /*
         * when `key` is an empty string retrieve the object meta-data of the bucket
         */
        if( "".equals(s3Path.getKey()) ) {
            ObjectMetadata meta = client.getObjectMetadata(s3Path.getBucket(), "");
            if( meta == null )
                throw new NoSuchFileException("s3://" + s3Path.getBucket());

            summary = new S3ObjectSummary();
            summary.setBucketName(s3Path.getBucket());
            summary.setETag(meta.getETag());
            summary.setKey(s3Path.getKey());
            summary.setLastModified(meta.getLastModified());
            summary.setSize(meta.getContentLength());
            // TODO summary.setOwner(?);
            // TODO summary.setStorageClass(?);
            return summary;
        }

        /*
         * Lookup for the object summary for the specified object key
         * by using a `listObjects` request
         */
        String marker = null;
        while( true ) {
            ListObjectsRequest request = new ListObjectsRequest();
            request.setBucketName(s3Path.getBucket());
            request.setPrefix(s3Path.getKey());
            request.setMaxKeys(250);
            if( marker != null )
                request.setMarker(marker);

            ObjectListing listing = client.listObjects(request);
            List<S3ObjectSummary> results = listing.getObjectSummaries();

            if (results.isEmpty()){
                break;
            }

            for( S3ObjectSummary item : results ) {
                if( matchName(s3Path.getKey(), item)) {
                    return item;
                }
            }

            if( listing.isTruncated() )
                marker = listing.getNextMarker();
            else
                break;
        }

        throw new NoSuchFileException("s3://" + s3Path.getBucket() + "/" + s3Path.getKey());
    }

    private boolean matchName(String fileName, S3ObjectSummary summary) {
        String foundKey = summary.getKey();

        // they are different names return false
        if( !foundKey.startsWith(fileName) ) {
            return false;
        }

        // when they are the same length, they are identical
        if( foundKey.length() == fileName.length() )
            return true;

        return foundKey.charAt(fileName.length()) == '/';
    }

    public ObjectMetadata getS3ObjectMetadata(S3Path s3Path) {
        S3Client client = s3Path.getFileSystem().getClient();
        try {
            return client.getObjectMetadata(s3Path.getBucket(), s3Path.getKey());
        }
        catch (AmazonS3Exception e){
            if (e.getStatusCode() != 404){
                throw e;
            }
            return null;
        }
    }

    /**
     * get S3Object represented by this S3Path try to access with or without end slash '/'
     * @param s3Path S3Path
     * @return S3Object or null if not exists
     */
    @Deprecated
    private S3Object getS3Object(S3Path s3Path){

        S3Client client = s3Path.getFileSystem()
                .getClient();

        S3Object object = getS3Object(s3Path.getBucket(), s3Path.getKey(), client);

        if (object != null) {
            return object;
        }
        else{
            return getS3Object(s3Path.getBucket(), s3Path.getKey() + "/", client);
        }
    }

    /**
     * get s3Object with S3Object#getObjectContent closed
     * @param bucket String bucket
     * @param key String key
     * @param client S3Client client
     * @return S3Object
     */
    private S3Object getS3Object(String bucket, String key, S3Client client){
        try {
            S3Object object = client .getObject(bucket, key);
            if (object.getObjectContent() != null){
                try {
                    object.getObjectContent().close();
                }
                catch (IOException e ) {
                    log.debug("Error while closing S3Object for bucket: `{}` and key: `{}` -- Cause: {}",bucket, key, e.getMessage());
                }
            }
            return object;
        }
        catch (AmazonS3Exception e){
            if (e.getStatusCode() != 404){
                throw e;
            }
            return null;
        }
    }
}
