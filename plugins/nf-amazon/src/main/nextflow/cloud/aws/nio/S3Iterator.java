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

package nextflow.cloud.aws.nio;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import software.amazon.awssdk.services.s3.model.*;
import com.google.common.base.Preconditions;

/**
 * S3 iterator over folders at first level.
 * Future versions of this class should be return the elements
 * in a incremental way when the #next() method is called.
 */
public class S3Iterator implements Iterator<Path> {

    private S3FileSystem s3FileSystem;
    private String bucket;
    private String key;

    private Iterator<S3Path> it;

    public S3Iterator(S3FileSystem s3FileSystem, String bucket, String key) {

        Preconditions.checkArgument(key != null && key.endsWith("/"), "key %s should be ended with slash '/'", key);

        this.bucket = bucket;
        // the only case i dont need the end slash is to list buckets content
        this.key = key.length() == 1 ? "" : key;
        this.s3FileSystem = s3FileSystem;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    @Override
    public S3Path next() {
        return getIterator().next();
    }

    @Override
    public boolean hasNext() {
        return getIterator().hasNext();
    }

    private Iterator<S3Path> getIterator() {
        if (it == null) {
            ListObjectsV2Request request = buildRequest();

            S3Client s3Client = s3FileSystem.getClient();

            // This automatically handles pagination
            it = s3Client.listObjectsV2Paginator(request).stream().flatMap(r -> parseObjectListing(r).stream()).iterator();
        }

        return it;
    }

    private ListObjectsV2Request buildRequest(){

        return ListObjectsV2Request.builder()
                .bucket(bucket)
                .prefix(key)
                .delimiter("/")
                .build();
    }

    /**
     * add to the listPath the elements at the same level that s3Path
     * @param current ListObjectsResponseto walk
     */
    private List<S3Path> parseObjectListing( ListObjectsV2Response current) {
        List<S3Path> listPath = new ArrayList<>();
        // add all the objects i.e. the files, except iterator key.
        // In V2, object listing is also returning the key of the request. Skip it from the iterator to avoid loops.
        for (final S3Object objectSummary : current.contents()) {
            final String key = objectSummary.key();
            if( this.key.equals(key)) continue;
            final S3Path path = new S3Path(s3FileSystem, "/" + bucket, key.split("/"));
            path.setObjectSummary(objectSummary);
            listPath.add(path);
        }

        // add all the common prefixes i.e. the directories, except iterator key
        for(final CommonPrefix prefix : current.commonPrefixes()) {
            if( prefix.prefix().equals("/") || this.key.equals(prefix.prefix())) continue;
            listPath.add(new S3Path(s3FileSystem, "/" + bucket, prefix.prefix()));
        }
        return listPath;
    }

    /**
     * The current #buildRequest() get all subdirectories and her content.
     * This method filter the keyChild and check if is a immediate
     * descendant of the keyParent parameter
     * @param keyParent String
     * @param keyChild String
     * @return String parsed
     *  or null when the keyChild and keyParent are the same and not have to be returned
     */
    @Deprecated
    private String getInmediateDescendent(String keyParent, String keyChild){

        keyParent = deleteExtraPath(keyParent);
        keyChild = deleteExtraPath(keyChild);

        final int parentLen = keyParent.length();
        final String childWithoutParent = deleteExtraPath(keyChild
                .substring(parentLen));

        String[] parts = childWithoutParent.split("/");

        if (parts.length > 0 && !parts[0].isEmpty()){
            return keyParent + "/" + parts[0];
        }
        else {
            return null;
        }

    }

    @Deprecated
    private String deleteExtraPath(String keyChild) {
        if (keyChild.startsWith("/")){
            keyChild = keyChild.substring(1);
        }
        if (keyChild.endsWith("/")){
            keyChild = keyChild.substring(0, keyChild.length() - 1);
        }
        return keyChild;
    }
}
