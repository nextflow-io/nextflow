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

package nextflow.cloud.aws.nio

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.file.ObjectMeta
import nextflow.file.ObjectStoreReader
import org.pf4j.Extension
import software.amazon.awssdk.services.s3.model.GetObjectRequest
import software.amazon.awssdk.services.s3.model.ListObjectsV2Request
import software.amazon.awssdk.services.s3.model.ListObjectsV2Response
import software.amazon.awssdk.services.s3.model.S3Object

/**
 * {@link ObjectStoreReader} for {@code s3://} paths.
 *
 * <ul>
 *   <li>a single no-delimiter {@code ListObjectsV2} yields every descendant's
 *       {@code (relpath, size, mtime)} in {@code O(N/1000)} calls, so a record-backed identity can
 *       guard/aggregate a directory without a per-object HEAD or a per-subfolder delimiter LIST;</li>
 *   <li>a single ranged {@code GetObject} (via the {@code Range} header) reads a byte span in one
 *       call, so the {@code sample} identity can hash a few small regions of a large object without
 *       downloading it in full.</li>
 * </ul>
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Extension
@CompileStatic
class S3ObjectStoreReader extends ObjectStoreReader {

    @Override
    boolean canHandle(String scheme) {
        return scheme == 's3'
    }

    @Override
    List<Map.Entry<String,ObjectMeta>> listWithMeta(Path prefix) {
        final s3 = (S3Path) prefix
        final client = s3.getFileSystem().getClient()
        final bucket = s3.getBucket()
        final base = normalizePrefix(s3.getKey())
        // NO delimiter -> a single flat/recursive listing of all descendants (paginated at 1000/page),
        // so the call count is O(N/1000) regardless of how many subfolders the tree has.
        final req = ListObjectsV2Request.builder().bucket(bucket).prefix(base).build()
        final out = new ArrayList<Map.Entry<String,ObjectMeta>>()
        for( ListObjectsV2Response page : client.listObjectsV2Paginator(req) ) {
            for( S3Object o : page.contents() ) {
                final rel = relativize(base, o.key())
                if( rel == null )
                    continue
                out.add(Map.entry(rel, new ObjectMeta(o.size(), o.lastModified().toEpochMilli())))
            }
        }
        return out
    }

    @Override
    byte[] readRange(Path path, long offset, int len) {
        final s3 = (S3Path) path
        final client = s3.getFileSystem().getClient()
        final req = GetObjectRequest.builder()
                .bucket(s3.getBucket())
                .key(s3.getKey())
                .range("bytes=${offset}-${offset + len - 1}".toString())
                .build()
        return client.getObjectRange(req)
    }
}
