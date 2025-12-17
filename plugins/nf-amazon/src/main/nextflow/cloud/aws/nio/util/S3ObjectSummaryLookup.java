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

package nextflow.cloud.aws.nio.util;

import java.io.IOException;
import java.nio.file.AccessDeniedException;
import java.nio.file.NoSuchFileException;
import java.util.List;

import nextflow.cloud.aws.nio.S3Client;
import nextflow.cloud.aws.nio.S3Path;
import software.amazon.awssdk.awscore.exception.AwsServiceException;
import software.amazon.awssdk.core.exception.SdkException;
import software.amazon.awssdk.services.s3.model.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class S3ObjectSummaryLookup {

    private static final Logger log = LoggerFactory.getLogger(S3Object.class);

    /**
     * Get the {@link software.amazon.awssdk.services.s3.model.S3Object} that represent this Path or its first child if the path does not exist
     * @param s3Path {@link S3Path}
     * @return {@link software.amazon.awssdk.services.s3.model.S3Object}
     * @throws java.nio.file.NoSuchFileException if not found the path and any child
     * @throws java.nio.file.AccessDeniedException if access is denied due to missing credentials or permissions
     * @throws java.io.IOException for other S3 access errors
     */
    public S3Object lookup(S3Path s3Path) throws IOException {

        /*
         * check is object summary has been cached
         */
        S3Object summary = s3Path.fetchObject();
        if( summary != null ) {
            return summary;
        }

        final S3Client client = s3Path.getFileSystem().getClient();
        final String path = s3Path.toS3ObjectId().toString();

        try {
            /*
             * when `key` is an empty string retrieve the object meta-data of the bucket
             */
            if( "".equals(s3Path.getKey()) ) {
                HeadObjectResponse meta = client.getObjectMetadata(s3Path.getBucket(), "");
                if( meta == null )
                    throw new NoSuchFileException("s3://" + s3Path.getBucket());

                summary = S3Object.builder()
                        .eTag(meta.eTag())
                        .key(s3Path.getKey())
                        .lastModified(meta.lastModified())
                        .size(meta.contentLength())
                        .build();

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
                ListObjectsRequest.Builder request = ListObjectsRequest.builder();
                request.bucket(s3Path.getBucket());
                request.prefix(s3Path.getKey());
                request.maxKeys(250);
                if( marker != null )
                    request.marker(marker);

                ListObjectsResponse listing = client.listObjects(request.build());
                List<S3Object> results = listing.contents();

                if (results.isEmpty()){
                    break;
                }

                for( S3Object item : results ) {
                    if( matchName(s3Path.getKey(), item)) {
                        return item;
                    }
                }

                if( listing.isTruncated() )
                    marker = listing.nextMarker();
                else
                    break;
            }

            throw new NoSuchFileException("s3://" + s3Path.getBucket() + "/" + s3Path.getKey());
        }
        catch (SdkException e) {
            throw translateException(e, path);
        }
    }

    /**
     * Translate AWS SDK exceptions to standard Java NIO exceptions with clear error messages.
     *
     * @param e the AWS SDK exception
     * @param path the S3 path being accessed (for error messages)
     * @return an appropriate IOException subclass
     */
    IOException translateException(SdkException e, String path) {
        // Check for service exceptions with HTTP status codes
        if (e instanceof AwsServiceException) {
            int statusCode = ((AwsServiceException) e).statusCode();
            if (statusCode == 403 || statusCode == 401) {
                return new AccessDeniedException(path, null,
                    "Access denied to S3 path - check AWS credentials and permissions");
            }
            if (statusCode == 404) {
                return new NoSuchFileException(path);
            }
        }

        // Check for client-side credential errors based on message patterns
        String message = e.getMessage();
        if (message != null && message.contains("Unable to load credentials")) {
            return new AccessDeniedException(path, null,
                "Cannot access S3 path - AWS credentials may not be configured");
        }

        // Default: wrap as IOException with original cause
        return new IOException("S3 error accessing path: " + path, e);
    }

    private boolean matchName(String fileName, S3Object summary) {
        String foundKey = summary.key();

        // they are different names return false
        if( !foundKey.startsWith(fileName) ) {
            return false;
        }

        // when they are the same length, they are identical
        if( foundKey.length() == fileName.length() )
            return true;

        return foundKey.charAt(fileName.length()) == '/';
    }

    public HeadObjectResponse getS3ObjectMetadata(S3Path s3Path) throws IOException {
        S3Client client = s3Path.getFileSystem().getClient();
        final String path = s3Path.toS3ObjectId().toString();
        try {
            return client.getObjectMetadata(s3Path.getBucket(), s3Path.getKey());
        }
        catch (SdkException e) {
            // Check if it's a 404 - return null for not found
            if (e instanceof AwsServiceException && ((AwsServiceException) e).statusCode() == 404) {
                return null;
            }
            throw translateException(e, path);
        }
    }
}
