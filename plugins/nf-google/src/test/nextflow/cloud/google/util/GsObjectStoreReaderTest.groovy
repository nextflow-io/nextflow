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
import java.time.Instant
import java.time.ZoneOffset

import com.google.api.gax.paging.Page
import com.google.cloud.ReadChannel
import com.google.cloud.storage.Blob
import com.google.cloud.storage.BlobId
import com.google.cloud.storage.Storage
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import nextflow.Global
import nextflow.Session
import spock.lang.Specification
import nextflow.file.ObjectMeta

class GsObjectStoreReaderTest extends Specification {

    // test subclass injecting a mocked Storage client (bypasses the real StorageOptions build)
    static class Local extends GsObjectStoreReader {
        Storage stub
        Local(Storage stub) { this.stub = stub }
        @Override protected Storage storage() { return stub }
    }

    private Blob blob(String name, long size, long mtimeMillis) {
        Stub(Blob) {
            getName() >> name
            getSize() >> size
            getUpdateTimeOffsetDateTime() >> Instant.ofEpochMilli(mtimeMillis).atOffset(ZoneOffset.UTC)
        }
    }

    def 'canHandle only gs'() {
        expect:
        new GsObjectStoreReader().canHandle('gs')
        !new GsObjectStoreReader().canHandle('s3')
    }

    def 'listWithMeta flat-lists descendants as (relpath, ObjectMeta), strips the prefix, skips placeholders'() {
        given:
        def storage = Mock(Storage)
        // a REAL CloudStoragePath (final classes can't be mocked); construction hits no network
        def path = CloudStorageFileSystem.forBucket('b').getPath('/dir')
        def page = Stub(Page) {
            iterateAll() >> [
                    blob('dir/', 0L, 1L),                 // placeholder -> skipped
                    blob('dir/out.bin', 10L, 1000L),
                    blob('dir/sub/a.bin', 20L, 2000L),
            ]
        }
        when:
        def result = new Local(storage).listWithMeta(path)
        then:
        1 * storage.list('b', _, _) >> page   // prefix + fields BlobListOptions
        result == [ Map.entry('out.bin', new ObjectMeta(10, 1000)), Map.entry('sub/a.bin', new ObjectMeta(20, 2000)) ]
    }

    def 'requests the NAME field in the projection (else Blob.getName() is null at runtime)'() {
        given:
        def storage = Mock(Storage)
        def path = CloudStorageFileSystem.forBucket('b').getPath('/dir')
        def page = Stub(Page) { iterateAll() >> [] }
        when:
        new Local(storage).listWithMeta(path)
        then: 'the second BlobListOption is the fields projection and it includes NAME (with SIZE + UPDATED)'
        1 * storage.list('b', _, { it == Storage.BlobListOption.fields(Storage.BlobField.NAME, Storage.BlobField.SIZE, Storage.BlobField.UPDATED) }) >> page
    }

    def 'a task-dir LIST uses a slash-terminated prefix, so `…/<hash>` never matches `…/<hash>-2`'() {
        given: 'the first attempt dir of a task -- its `-2` sibling shares every byte of the key'
        def storage = Mock(Storage)
        def path = CloudStorageFileSystem.forBucket('b').getPath('/cache/work/ab/cdef0123456789abcdef0123456789')
        def page = Stub(Page) { iterateAll() >> [] }
        when:
        new Local(storage).listWithMeta(path)
        then: 'the prefix option ends with a slash -- without it the sibling attempt`s objects would be listed as members'
        1 * storage.list('b', { it == Storage.BlobListOption.prefix('cache/work/ab/cdef0123456789abcdef0123456789/') }, _) >> page
    }

    def 'readRange seeks to the offset and reads len bytes from the blob'() {
        given:
        def storage = Mock(Storage)
        def path = CloudStorageFileSystem.forBucket('b').getPath('/dir/file.bin')
        def data = [1,2,3,4,5] as byte[]
        def reader = Mock(ReadChannel)
        when:
        def result = new Local(storage).readRange(path, 100, 5)
        then:
        1 * storage.reader({ BlobId id -> id.getBucket() == 'b' && id.getName() == 'dir/file.bin' }) >> reader
        1 * reader.seek(100)
        and: 'the read is BOUNDED both ways -- otherwise the SDK fetches and allocates its 2 MiB'
        and: 'default chunk for a 5-byte window, and `limit` takes an absolute end position'
        1 * reader.limit(105) >> reader
        1 * reader.setChunkSize(5)
        and:
        1 * reader.read(_ as ByteBuffer) >> { ByteBuffer buf -> buf.put(data); return data.length }
        1 * reader.close()
        result == data
    }

    def 'the session config is re-read per call, so a second Session does not get the first one\'s'() {
        given: 'one extension instance, as SingletonExtensionFactory creates it: it lives for the JVM'
        def provider = new GsObjectStoreReader()

        when: 'the first session bills a requester-pays bucket to project-a'
        Global.session = Stub(Session) {
            getConfig() >> [google: [project: 'project-a', enableRequesterPaysBuckets: true]]
        }
        then:
        provider.userProject() == 'project-a'

        when: 'a second session replaces it -- tests, or embedded use'
        Global.session = Stub(Session) {
            getConfig() >> [google: [project: 'project-b', enableRequesterPaysBuckets: true]]
        }
        then: 'the new config wins; an instance memo here would have pinned project-a for the JVM'
        provider.userProject() == 'project-b'

        when: 'requester-pays is off, nothing is billed'
        Global.session = Stub(Session) { getConfig() >> [google: [project: 'project-c']] }
        then:
        provider.userProject() == null

        cleanup:
        Global.session = null
    }
}
