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

package nextflow.cloud.google.util

import java.nio.ByteBuffer
import java.nio.file.Path

import com.google.cloud.ReadChannel
import com.google.cloud.storage.Blob
import com.google.cloud.storage.BlobId
import com.google.cloud.storage.Storage
import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.CompileStatic
import nextflow.file.ObjectMeta
import nextflow.file.ObjectStoreReader
import org.pf4j.Extension

/**
 * {@link ObjectStoreReader} for {@code gs://} paths.
 *
 * <ul>
 *   <li>a single prefix {@code Storage.list} (no delimiter, requesting only the NAME, SIZE and
 *       UPDATED fields) yields every descendant's {@code (relpath, size, mtime)} in
 *       {@code O(N/1000)} calls, so a record-backed identity can guard/aggregate a directory without
 *       a per-object {@code get} — this is what makes GCS directories (and GCS task outputs,
 *       whose record is built from the same listing) content-addressable;</li>
 *   <li>a single seeked {@code Storage.reader} read fetches a byte span in one call, so the
 *       {@code sample} identity can hash a few small regions of a large object without downloading
 *       it in full.</li>
 * </ul>
 *
 * <p>Google's {@code CloudStorageFileSystemProvider} cannot be extended, so both halves use the GCS
 * {@code Storage} client directly rather than going through NIO.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Extension
@CompileStatic
class GsObjectStoreReader extends ObjectStoreReader {

    @Override
    boolean canHandle(String scheme) {
        return scheme == 'gs'
    }

    /**
     * The GCS client, carrying the configured credentials, project id, transport timeouts and retry
     * policy — see {@link GsStorageOptions}. One per instance, and this class is one extension, so
     * merging the listing and range halves also merged their two independent caches into one client.
     */
    protected Storage storage() {
        // NOT @Memoized: this extension is instantiated once per JVM (SingletonExtensionFactory) and
        // lives in the static providers list, so an instance memo would outlive the session and hand
        // the NEXT one the first session's credentials, project and requester-pays setting -- the
        // very thing GsStorageOptions.sessionOpts warns about. `sharedClientFor` is memoized keyed by
        // the options, so the client is still built once per config; this call is a map lookup.
        return GsStorageOptions.sharedClientFor(GsStorageOptions.sessionOpts())
    }

    /** The project to bill on a requester-pays bucket, or {@code null}. Not memoized -- see {@link #storage}. */
    protected String userProject() {
        return GsStorageOptions.userProject(GsStorageOptions.sessionOpts())
    }

    /** The object key of a GCS path: its string form without the leading {@code '/'}. */
    private static String keyOf(CloudStoragePath gs) {
        final str = gs.toString()
        return str.startsWith('/') ? str.substring(1) : str
    }

    @Override
    List<Map.Entry<String,ObjectMeta>> listWithMeta(Path prefix) {
        final gs = (CloudStoragePath) prefix
        final bucket = gs.getFileSystem().bucket()
        final base = normalizePrefix(keyOf(gs))
        // NO delimiter -> a single flat/recursive listing of all descendants (paginated at 1000/page);
        // request only the fields we need (name + size + last-updated) to keep the payload small.
        // NAME is mandatory: with a fields() projection, getName() is null unless explicitly requested.
        final opts = new ArrayList<Storage.BlobListOption>()
        opts.add(Storage.BlobListOption.prefix(base))
        opts.add(Storage.BlobListOption.fields(Storage.BlobField.NAME, Storage.BlobField.SIZE, Storage.BlobField.UPDATED))
        final billTo = userProject()
        if( billTo )
            opts.add(Storage.BlobListOption.userProject(billTo))
        final page = storage().list(bucket, opts as Storage.BlobListOption[])
        final out = new ArrayList<Map.Entry<String,ObjectMeta>>()
        for( Blob b : page.iterateAll() ) {
            final rel = relativize(base, b.getName())
            if( rel == null )
                continue
            final mtime = b.getUpdateTimeOffsetDateTime().toInstant().toEpochMilli()
            out.add(Map.entry(rel, new ObjectMeta(b.getSize(), mtime)))
        }
        return out
    }

    @Override
    protected byte[] readRange0(Path path, long offset, int len) {
        final gs = (CloudStoragePath) path
        final bucket = gs.getFileSystem().bucket()
        // not try-with-resources: `limit` may return a DIFFERENT channel per its contract, and the
        // resource variable would be implicitly final. Close whichever instance we ended up reading
        // from -- in the current SDK that is the same object, but the contract does not promise it.
        final billTo = userProject()
        ReadChannel reader = billTo
                ? storage().reader(BlobId.of(bucket, keyOf(gs)), Storage.BlobSourceOption.userProject(billTo))
                : storage().reader(BlobId.of(bucket, keyOf(gs)))
        try {
            reader.seek(offset)
            // BOUND THE READ, twice over. Without this the channel keeps the SDK's default 2 MiB
            // chunk: it allocates that much and fetches it from the server, so each 16 KiB sample
            // window pulled 2 MiB -- 6 MiB per sampled file instead of 48 KiB, three times over.
            // That defeats the whole premise of this identity, which is that identifying a 200 GB
            // object costs the same as identifying a 2 MB one.
            //
            //  - `limit` bounds the wire request. It takes an ABSOLUTE end position, not a length
            //    (the range is [seek, limit)), and the contract says to use the channel it returns.
            //  - `setChunkSize` bounds the buffer the channel allocates.
            // Either one alone leaves half the waste: limit still allocates 2 MiB locally,
            // setChunkSize still leaves the request open to the end of the object.
            reader = reader.limit(offset + (long) len)
            reader.setChunkSize(len)
            final buf = ByteBuffer.allocate(len)
            while( buf.hasRemaining() && reader.read(buf) > 0 ) { }
            return Arrays.copyOf(buf.array(), buf.position())
        }
        finally {
            reader.close()
        }
    }
}
