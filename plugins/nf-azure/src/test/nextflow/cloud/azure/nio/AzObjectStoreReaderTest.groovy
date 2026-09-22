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

package nextflow.cloud.azure.nio

import java.time.Instant
import java.time.ZoneOffset

import com.azure.core.http.HttpMethod
import com.azure.core.http.HttpRequest
import com.azure.core.util.BinaryData
import com.azure.storage.blob.models.BlobDownloadContentAsyncResponse
import com.azure.storage.blob.models.BlobDownloadContentResponse
import com.azure.storage.blob.models.BlobItem
import com.azure.storage.blob.models.BlobItemProperties
import com.azure.storage.blob.models.BlobRange
import spock.lang.Specification
import nextflow.file.ObjectMeta

class AzObjectStoreReaderTest extends Specification {

    /** Test subclass injecting both seams directly (bypasses the real clients / network). */
    static class Local extends AzObjectStoreReader {
        List<BlobItem> items
        String capturedPrefix
        BlobDownloadContentResponse resp
        BlobRange capturedRange
        @Override protected Iterable<BlobItem> listBlobs(AzPath dir, String base) { capturedPrefix = base; return items }
        @Override protected BlobDownloadContentResponse downloadContent(AzPath az, BlobRange range) { capturedRange = range; return resp }
    }

    private BlobItem blob(String name, long size, long mtimeMillis) {
        new BlobItem()
                .setName(name)
                .setProperties(new BlobItemProperties()
                        .setContentLength(size)
                        .setLastModified(Instant.ofEpochMilli(mtimeMillis).atOffset(ZoneOffset.UTC)))
    }

    def 'canHandle only az'() {
        expect:
        new AzObjectStoreReader().canHandle('az')
        !new AzObjectStoreReader().canHandle('s3')
    }

    def 'listWithMeta flat-lists descendants as (relpath, ObjectMeta), strips the prefix, skips placeholders'() {
        given:
        def path = Mock(AzPath) { blobName() >> 'dir' }
        def provider = new Local(items: [
                blob('dir/', 0L, 1L),                 // placeholder -> skipped
                blob('dir/out.bin', 10L, 1000L),
                blob('dir/sub/a.bin', 20L, 2000L),
        ])
        when:
        def result = provider.listWithMeta(path)
        then: 'the prefix is normalized with a trailing slash and members are relativized against it'
        provider.capturedPrefix == 'dir/'
        result == [ Map.entry('out.bin', new ObjectMeta(10, 1000)), Map.entry('sub/a.bin', new ObjectMeta(20, 2000)) ]
    }

    def 'listWithMeta skips the Azure empty-dir marker, which is a placeholder and not a member'() {
        given: 'a published tree as NIO leaves it: copying into a subfolder materialises it with a marker'
        def path = Mock(AzPath) { blobName() >> 'pub/outdir' }
        def provider = new Local(items: [
                blob('pub/outdir/a.bin', 10L, 1000L),
                blob('pub/outdir/sub/' + AzFileSystem.EMPTY_DIR_MARKER, 0L, 1500L),
                blob('pub/outdir/sub/b.bin', 20L, 2000L),
                blob('pub/outdir/' + AzFileSystem.EMPTY_DIR_MARKER, 0L, 900L),
        ])

        when:
        def result = provider.listWithMeta(path)

        then: 'only real content is a member -- S3 drops its "<key>/" equivalent and GCS creates none,'
        and: 'so without this the same tree would enumerate differently per cloud and no recorded'
        and: 'directory identity could ever match the listing it is guarded against'
        result == [ Map.entry('a.bin', new ObjectMeta(10, 1000)), Map.entry('sub/b.bin', new ObjectMeta(20, 2000)) ]
    }

    def 'listWithMeta skips an ADLS Gen2 directory blob, which is a real blob but not a member'() {
        given: 'a hierarchical-namespace account, where a folder IS a blob: hdi_isfolder=true, size 0'
        def path = Mock(AzPath) { blobName() >> 'dir' }
        def provider = new Local(items: [
                blob('dir/sub', 0L, 1500L).setMetadata([hdi_isfolder: 'true']),
                blob('dir/out.bin', 10L, 1000L),
                blob('dir/sub/a.bin', 20L, 2000L),
        ])

        when:
        def result = provider.listWithMeta(path)

        then: 'the folder blob is not a member -- it has no content, and its mtime moves whenever a'
        and: 'child is added, so recording it would invalidate the directory guard on every write'
        result == [ Map.entry('out.bin', new ObjectMeta(10, 1000)), Map.entry('sub/a.bin', new ObjectMeta(20, 2000)) ]
    }

    def 'a task-dir LIST uses a slash-terminated prefix, so `…/<hash>` never matches `…/<hash>-2`'() {
        given: 'the first attempt dir of a task -- its `-2` sibling shares every byte of the blob name'
        def path = Mock(AzPath) { blobName() >> 'cache/work/ab/cdef0123456789abcdef0123456789' }
        def provider = new Local(items: [ blob('cache/work/ab/cdef0123456789abcdef0123456789/out.bin', 10L, 1000L) ])
        when:
        def result = provider.listWithMeta(path)
        then: 'the prefix ends with a slash -- without it the sibling attempt`s objects would be listed as members'
        provider.capturedPrefix == 'cache/work/ab/cdef0123456789abcdef0123456789/'
        result == [ Map.entry('out.bin', new ObjectMeta(10, 1000)) ]
    }

    def 'readRange requests the (offset, len) range and passes through the returned bytes'() {
        given:
        def path = Mock(AzPath)
        def bytes = [1,2,3] as byte[]
        // BlobDownloadContentResponse is final (can't be mocked/stubbed) -> build a real instance
        def asyncResp = new BlobDownloadContentAsyncResponse(
                new HttpRequest(HttpMethod.GET, "http://example.com"), 200, null, BinaryData.fromBytes(bytes), null)
        def provider = new Local(resp: new BlobDownloadContentResponse(asyncResp))
        when:
        def result = provider.readRange(path, 100, 5)
        then:
        provider.capturedRange.getOffset() == 100L
        provider.capturedRange.getCount() == 5L
        result == bytes
    }
}
