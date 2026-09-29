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

import java.time.Instant

import software.amazon.awssdk.services.s3.model.ListObjectsV2Response
import software.amazon.awssdk.services.s3.model.S3Object
import software.amazon.awssdk.services.s3.paginators.ListObjectsV2Iterable
import spock.lang.Specification
import nextflow.file.ObjectMeta

class S3ObjectStoreReaderTest extends Specification {

    def 'canHandle only s3'() {
        expect:
        new S3ObjectStoreReader().canHandle('s3')
        !new S3ObjectStoreReader().canHandle('gs')
    }

    def 'listWithMeta flat-lists descendants as (relpath, ObjectMeta), strips the prefix, skips placeholders'() {
        given:
        def client = Mock(S3Client)
        def fs = Mock(S3FileSystem) { getClient() >> client }
        def path = Mock(S3Path) { getFileSystem() >> fs; getBucket() >> 'b'; getKey() >> 'dir' }
        def page = ListObjectsV2Response.builder().contents(
                S3Object.builder().key('dir/').size(0L).lastModified(Instant.ofEpochMilli(1)).build(),        // placeholder -> skipped
                S3Object.builder().key('dir/out.bin').size(10L).lastModified(Instant.ofEpochMilli(1000)).build(),
                S3Object.builder().key('dir/sub/a.bin').size(20L).lastModified(Instant.ofEpochMilli(2000)).build(),
        ).build()
        when:
        def result = new S3ObjectStoreReader().listWithMeta(path)
        then:
        1 * client.listObjectsV2Paginator({ it.prefix() == 'dir/' }) >> Stub(ListObjectsV2Iterable) { iterator() >> [page].iterator() }
        result == [ Map.entry('out.bin', new ObjectMeta(10, 1000)), Map.entry('sub/a.bin', new ObjectMeta(20, 2000)) ]
    }

    def 'a task-dir LIST uses a slash-terminated prefix, so `…/<hash>` never matches `…/<hash>-2`'() {
        given: 'the first attempt dir of a task -- its `-2` sibling shares every byte of the key'
        def client = Mock(S3Client)
        def fs = Mock(S3FileSystem) { getClient() >> client }
        def path = Mock(S3Path) { getFileSystem() >> fs; getBucket() >> 'b'; getKey() >> 'cache/work/ab/cdef0123456789abcdef0123456789' }
        def page = ListObjectsV2Response.builder().contents(
                S3Object.builder().key('cache/work/ab/cdef0123456789abcdef0123456789/out.bin').size(10L).lastModified(Instant.ofEpochMilli(1000)).build(),
        ).build()
        when:
        def result = new S3ObjectStoreReader().listWithMeta(path)
        then: 'the prefix ends with a slash -- without it the sibling attempt`s objects would be listed as members'
        1 * client.listObjectsV2Paginator({ it.prefix() == 'cache/work/ab/cdef0123456789abcdef0123456789/' }) >> Stub(ListObjectsV2Iterable) { iterator() >> [page].iterator() }
        result == [ Map.entry('out.bin', new ObjectMeta(10, 1000)) ]
    }

    def 'readRange issues a single ranged GetObject and passes through the returned bytes'() {
        given:
        def client = Mock(S3Client)
        def fs = Mock(S3FileSystem) { getClient() >> client }
        def path = Mock(S3Path) { getFileSystem() >> fs; getBucket() >> 'b'; getKey() >> 'dir/file.bin' }
        def bytes = [1,2,3] as byte[]
        when:
        def result = new S3ObjectStoreReader().readRange(path, 100, 16384)
        then:
        1 * client.getObjectRange({ it.bucket() == 'b' && it.key() == 'dir/file.bin' && it.range() == 'bytes=100-16483' }) >> bytes
        result.is(bytes)
    }
}
