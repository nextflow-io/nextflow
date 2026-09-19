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

import java.nio.file.Path

import com.azure.storage.blob.models.BlobDownloadContentResponse
import com.azure.storage.blob.models.BlobItem
import com.azure.storage.blob.models.BlobRange
import com.azure.storage.blob.models.DownloadRetryOptions
import com.azure.storage.blob.models.ListBlobsOptions
import groovy.transform.CompileStatic
import nextflow.file.ObjectMeta
import nextflow.file.ObjectStoreReader
import org.pf4j.Extension

/**
 * {@link ObjectStoreReader} for {@code az://} paths.
 *
 * <ul>
 *   <li>a single no-delimiter {@code listBlobs} yields every descendant's
 *       {@code (relpath, size, mtime)} in {@code O(N/1000)} calls, so a record-backed identity can
 *       guard/aggregate an Azure directory without a per-blob HEAD;</li>
 *   <li>a single ranged blob download (via {@link BlobRange}) reads a byte span in one call, so the
 *       {@code sample} identity can hash a few small regions of a large blob without downloading it
 *       in full.</li>
 * </ul>
 *
 * <p>Azure already works without the listing half — {@code az} is a "cheaply listable" scheme, so the
 * NIO walk fallback serves member attrs from the bulk listing — but that walk issues one
 * <b>delimited</b> ({@code listBlobsByHierarchy}) LIST per subfolder. This collapses a whole subtree
 * into one flat listing, matching {@code S3ObjectStoreReader}/{@code GsObjectStoreReader}. Both paths
 * read the same {@code getContentLength()}/{@code getLastModified()} blob properties, so the recorded
 * {@code (size, mtime)} guard is identical regardless of which one produced it.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Extension
@CompileStatic
class AzObjectStoreReader extends ObjectStoreReader {

    @Override
    boolean canHandle(String scheme) {
        return scheme == 'az'
    }

    @Override
    List<Map.Entry<String,ObjectMeta>> listWithMeta(Path prefix) {
        final az = (AzPath) prefix
        final base = normalizePrefix(az.blobName())
        final out = new ArrayList<Map.Entry<String,ObjectMeta>>()
        for( BlobItem item : listBlobs(az, base) ) {
            final rel = relativize(base, item.getName())
            if( rel == null )
                continue
            // Azure materialises a folder with a placeholder object. It is not a member: S3's
            // equivalent (a key ending in '/') is already dropped by relativize() and GCS creates
            // none at all, so without this the SAME tree would enumerate differently per cloud --
            // and a directory identity recorded from one listing would never match another.
            if( rel == AzFileSystem.EMPTY_DIR_MARKER || rel.endsWith('/' + AzFileSystem.EMPTY_DIR_MARKER) )
                continue
            final size = item.getProperties().getContentLength()
            final mtime = item.getProperties().getLastModified().toInstant().toEpochMilli()
            out.add(Map.entry(rel, new ObjectMeta(size, mtime)))
        }
        return out
    }

    @Override
    byte[] readRange(Path path, long offset, int len) {
        final az = (AzPath) path
        final range = new BlobRange(offset, len as Long)
        final resp = downloadContent(az, range)
        return resp.getValue().toBytes()
    }

    /**
     * One flat (NO delimiter) listing of every descendant under {@code base}, so the call count is
     * {@code O(N/1000)} regardless of subfolder nesting. Declared {@code Iterable} (not the SDK's
     * {@code PagedIterable}) so a test can inject a plain list.
     */
    protected Iterable<BlobItem> listBlobs(AzPath dir, String base) {
        return dir.containerClient().listBlobs(new ListBlobsOptions().setPrefix(base), null)
    }

    /**
     * Issue the ranged download. Declared as a protected method (rather than inlined) so a test can
     * inject a stub blob client response without a real Azure connection.
     */
    protected BlobDownloadContentResponse downloadContent(AzPath az, BlobRange range) {
        return az.blobClient().downloadContentWithResponse(new DownloadRetryOptions(), null, range, false, null, null)
    }
}
