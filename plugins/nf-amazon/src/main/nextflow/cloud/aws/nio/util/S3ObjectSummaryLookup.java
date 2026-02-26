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
     * Lookup for the S3 object matching the specified path using at most two bounded
     * {@code listObjects} calls (replaces the previous unbounded pagination loop).
     *
     * @param s3Path the S3 path to look up
     * @param client the S3 client
     * @return the matching {@link S3Object}, or {@code null} if not found
     */
    private S3Object getS3Object(S3Path s3Path, S3Client client) {

        // Call 1: list up to 2 objects whose key starts with the target key.
        //
        // Why maxKeys(2) instead of paginating all results?
        // The previous implementation used an unbounded while(true) loop fetching 250 keys
        // per page. On prefixes with millions of objects this caused excessive S3 LIST API
        // calls, high latency, and potential timeouts. Two results are enough to cover
        // the common cases:
        //   - Exact file match: the key itself exists as an object (e.g. "data.txt")
        //   - Directory match: a child object (e.g. "data/file1") appears within the
        //     first 2 lexicographic results
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

        // Call 2 (fallback): list 1 object with prefix "key/" to detect directories
        // that Call 1 missed.
        //
        // Why can Call 1 miss a directory?
        // S3 lists keys in lexicographic (UTF-8 byte) order, and several common characters
        // sort *before* '/' (0x2F) — notably '-' (0x2D) and '.' (0x2E).
        //
        // Example: given keys "a-a/file-3", "a.txt", and "a/file-1", S3 returns them as:
        //   a-a/file-3   ← '-' (0x2D) < '/' (0x2F)
        //   a.txt         ← '.' (0x2E) < '/' (0x2F)
        //   a/file-1      ← '/' (0x2F) — the actual directory child
        //
        // With maxKeys(2), Call 1 only sees "a-a/file-3" and "a.txt" — neither matches
        // key "a" via matchName(). The directory child "a/file-1" is pushed beyond the
        // result window by sibling keys that sort earlier.
        //
        // By searching with prefix "a/" directly, we skip all those siblings and find
        // "a/file-1", confirming that "a" is a directory.
        request = ListObjectsRequest.builder()
                .bucket(s3Path.getBucket())
                .prefix(s3Path.getKey()+'/')
                .maxKeys(1)
                .build();

        listing = client.listObjects(request);
        results = listing.contents();
        for( S3Object item : results ) {
            if( matchName(s3Path.getKey(), item)) {
                return item;
            }
        }
        return null;
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
