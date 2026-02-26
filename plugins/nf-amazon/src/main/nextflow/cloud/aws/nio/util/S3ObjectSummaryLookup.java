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

import java.nio.file.NoSuchFileException;
import java.util.List;

import nextflow.cloud.aws.nio.S3Client;
import software.amazon.awssdk.services.s3.model.*;
import nextflow.cloud.aws.nio.S3Path;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class S3ObjectSummaryLookup {

    private static final Logger log = LoggerFactory.getLogger(S3Object.class);

    /**
     * Get the {@link software.amazon.awssdk.services.s3.model.S3Object} that represent this Path or its first child if the path does not exist
     * @param s3Path {@link S3Path}
     * @return {@link software.amazon.awssdk.services.s3.model.S3Object}
     * @throws java.nio.file.NoSuchFileException if not found the path and any child
     */
    public S3Object lookup(S3Path s3Path) throws NoSuchFileException {

        /*
         * check is object summary has been cached
         */
        S3Object summary = s3Path.fetchObject();
        if( summary != null ) {
            return summary;
        }

        final S3Client client = s3Path.getFileSystem().getClient();

        /*
         * when `key` is an empty string retrieve the object meta-data of the bucket
         */
        if( "".equals(s3Path.getKey()) ) {
            HeadBucketResponse meta = client.getBucketMetadata(s3Path.getBucket());
            if( meta == null )
                throw new NoSuchFileException("s3://" + s3Path.getBucket());

            summary = S3Object.builder()
                    .key(s3Path.getKey())
                    .build();

            // TODO summary.setOwner(?);
            // TODO summary.setStorageClass(?);
            return summary;
        }

        S3Object item = getS3Object(s3Path, client);
        if( item != null )
            return item;

        throw new NoSuchFileException("s3://" + s3Path.getBucket() + "/" + s3Path.getKey());
    }

    /**
     * Lookup for the S3 object matching the specified path using at one head request
     * and one bounded list-objects call.
     *
     * @param s3Path the S3 path to look up
     * @param client the S3 client
     * @return the matching {@link S3Object}, or {@code null} if not found
     */
    private S3Object getS3Object(S3Path s3Path, S3Client client) {
        /*
         * First: try HEAD request for exact object (fast, works on all bucket types including Express)
         */
        try {
            HeadObjectResponse metadata = client.getObjectMetadata(s3Path.getBucket(), s3Path.getKey());
            if( metadata != null ) {
                return S3Object.builder()
                        .key(s3Path.getKey())
                        .size(metadata.contentLength())
                        .lastModified(metadata.lastModified())
                        .eTag(metadata.eTag())
                        .storageClass(metadata.storageClassAsString())
                        .build();
            }
        }
        catch (S3Exception e) {
            if( e.statusCode() != 404 ) {
                throw e;
            }
            // 404 = object doesn't exist as a file, fall through to directory check
        }

        /*
         * Second: check if it's a "directory" by listing with trailing slash.
         * S3 Express One Zone buckets require prefixes to end with a delimiter.
         */
        String prefix = s3Path.getKey();
        if( !prefix.endsWith("/") ) {
            prefix = prefix + "/";
        }

        ListObjectsV2Request request = ListObjectsV2Request.builder()
                .bucket(s3Path.getBucket())
                .prefix(prefix)
                .maxKeys(1)  // only need to find one child to confirm directory exists
                .build();

        for( S3Object item : client.listObjectsV2Paginator(request).contents() ) {
            // Found a child, so this "directory" exists - return a synthetic S3Object for it
            return S3Object.builder()
                    .key(s3Path.getKey())
                    .build();
        }

        return null;
    }
}
