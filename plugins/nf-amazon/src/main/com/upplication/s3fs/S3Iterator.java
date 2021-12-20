/*
 * Copyright 2020, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 Javier Arnáiz @arnaix
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package com.upplication.s3fs;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import com.amazonaws.services.s3.model.ListObjectsRequest;
import com.amazonaws.services.s3.model.ObjectListing;
import com.amazonaws.services.s3.model.S3ObjectSummary;
import com.google.common.base.Preconditions;

/**
 * S3 iterator over folders at first level.
 * Future verions of this class should be return the elements
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
            List<S3Path> listPath = new ArrayList<>();

            // iterator over this list
            ObjectListing current = s3FileSystem.getClient().listObjects(buildRequest());

            while (current.isTruncated()) {
                // parse the elements
                parseObjectListing(listPath, current);
                // continue
                current = s3FileSystem.getClient().listNextBatchOfObjects(current);
            }

            parseObjectListing(listPath, current);

            it = listPath.iterator();
        }

        return it;
    }

    private ListObjectsRequest buildRequest(){

        ListObjectsRequest request = new ListObjectsRequest();
        request.setBucketName(bucket);
        request.setPrefix(key);
        request.setMarker(key);
        request.setDelimiter("/");
        return request;
    }

    /**
     * add to the listPath the elements at the same level that s3Path
     * @param listPath List not null list to add
     * @param current ObjectListing to walk
     */
    private void parseObjectListing(List<S3Path> listPath, ObjectListing current) {

        // add all the objects i.e. the files
        for (final S3ObjectSummary objectSummary : current.getObjectSummaries()) {
            final String key = objectSummary.getKey();
            final S3Path path = new S3Path(s3FileSystem, "/" + bucket, key.split("/"));
            path.setObjectSummary(objectSummary);
            listPath.add(path);
        }

        // add all the common prefixes i.e. the directories
        for(final String dir : current.getCommonPrefixes()) {
            if( dir.equals("/") ) continue;
            listPath.add(new S3Path(s3FileSystem, "/" + bucket, dir));
        }

    }

    /**
     * The current #buildRequest() get all subdirectories and her content.
     * This method filter the keyChild and check if is a inmediate
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
