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

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class S3UploadHelper {

    private static final long _1_KiB = 1024;
    private static final long _1_MiB = _1_KiB * _1_KiB;
    private static final long _1_GiB = _1_KiB * _1_KiB * _1_KiB;
    private static final long _1_TiB = _1_KiB * _1_KiB * _1_KiB * _1_KiB;

    /**
     * AWS S3 max part size
     * https://docs.aws.amazon.com/AmazonS3/latest/userguide/qfacts.html
     */
    public static final long MIN_PART_SIZE = 5 * _1_MiB;

    /**
     * AWS S3 min part size
     * https://docs.aws.amazon.com/AmazonS3/latest/userguide/qfacts.html
     */
    public static final long MAX_PART_SIZE = 5 * _1_GiB;

    /**
     * AWS S3 max object size
     * https://docs.aws.amazon.com/AmazonS3/latest/userguide/qfacts.html
     */
    public static final long MAX_OBJECT_SIZE = 5 * _1_TiB;

    /**
     * AWS S3 max parts in multi-part upload and copy request
     */
    public static final int MAX_PARTS_COUNT = 10_000;

    static public long computePartSize( long objectSize, long chunkSize ) {
        if( objectSize<0 ) throw new IllegalArgumentException("Argument 'objectSize' cannot be less than zero");
        if( chunkSize<MIN_PART_SIZE ) throw new IllegalArgumentException("Argument 'chunkSize' cannot be less than " + MIN_PART_SIZE);
        // Multipart upload and copy allows max 10_000 parts
        // each part can be up to 5 GB
        // Max file size is 5 TB
        // See https://docs.aws.amazon.com/AmazonS3/latest/userguide/qfacts.html
        long numOfParts = objectSize / chunkSize;
        if( numOfParts > MAX_PARTS_COUNT) {
            final long x = ceilDiv(objectSize, MAX_PARTS_COUNT);
            return ceilDiv(x, 10* _1_MiB) *10* _1_MiB;
        }
        return chunkSize;
    }


    private static long ceilDiv(long x, long y){
        return -Math.floorDiv(-x,y);
    }

    private static long ceilDiv(long x, int y){
        return -Math.floorDiv(-x,y);
    }

    static public void checkPartSize(long partSize) {
        if( partSize<MIN_PART_SIZE ) {
            String msg = String.format("The minimum part size for S3 multipart copy and upload operation cannot be less than 5 MiB -- offending value: %d", partSize);
            throw new IllegalArgumentException(msg);
        }

        if( partSize>MAX_PART_SIZE ) {
            String msg = String.format("The minimum part size for S3 multipart copy and upload operation cannot be less than 5 GiB -- offending value: %d", partSize);
            throw new IllegalArgumentException(msg);
        }
    }

    static public void checkPartIndex(int i, String path, long fileSize, long chunkSize) {
        if( i < 1 ) {
            String msg = String.format("S3 multipart copy request index cannot less than 1 -- offending value: %d; file: '%s'; size: %d; part-size: %d", i, path, fileSize, chunkSize);
            throw new IllegalArgumentException(msg);
        }
        if( i > MAX_PARTS_COUNT) {
            String msg = String.format("S3 multipart copy request exceed the number of max allowed parts -- offending value: %d; file: '%s'; size: %d; part-size: %d", i, path, fileSize, chunkSize);
            throw new IllegalArgumentException(msg);
        }
    }

}
