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

        /*
         * Lookup for the object summary for the specified object key
         * by using a `listObjects` request.
         *
         * Only 2 results are needed: the exact key match or its first
         * child (key + "/"). Fetching more would cause unbounded pagination
         * on prefixes with millions of objects.
         */
        ListObjectsRequest request = ListObjectsRequest.builder()
                .bucket(s3Path.getBucket())
                .prefix(s3Path.getKey())
                .maxKeys(2)
                .build();

        ListObjectsResponse listing = client.listObjects(request);
        List<S3Object> results = listing.contents();

        for( S3Object item : results ) {
            if( matchName(s3Path.getKey(), item)) {
                return item;
            }
        }

        throw new NoSuchFileException("s3://" + s3Path.getBucket() + "/" + s3Path.getKey());
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

    public HeadObjectResponse getS3ObjectMetadata(S3Path s3Path) {
        S3Client client = s3Path.getFileSystem().getClient();
        try {
            return client.getObjectMetadata(s3Path.getBucket(), s3Path.getKey());
        }
        catch (S3Exception e){
            if (e.statusCode() != 404){
                throw e;
            }
            return null;
        }
    }
}
