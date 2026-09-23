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

package nextflow.file

import java.nio.file.Path
import java.nio.file.Paths

import spock.lang.Specification

class ObjectStoreReaderTest extends Specification {

    static class FakeReader extends ObjectStoreReader {
        boolean canHandle(String scheme) { scheme == 'file' }
        List<Map.Entry<String,ObjectMeta>> listWithMeta(Path prefix) { return [Map.entry('a.bin', new ObjectMeta(10, 1))] }
        protected byte[] readRange0(Path path, long offset, int len) { return [1,2,3] as byte[] }
    }

    def 'lookup resolves a provider by scheme'() {
        given:
        def provider = new FakeReader()
        ObjectStoreReader.setProviders([provider] as List<ObjectStoreReader>)
        expect:
        ObjectStoreReader.lookup(Paths.get('/tmp/x')).is(provider)   // 'file' scheme
        cleanup:
        ObjectStoreReader.setProviders(null)
    }

    def 'lookup returns null when no provider handles the scheme'() {
        given: 'a read capability always has a documented fallback, so this must not throw'
        ObjectStoreReader.setProviders([] as List<ObjectStoreReader>)
        expect:
        ObjectStoreReader.lookup(Paths.get('/tmp/x')) == null
        cleanup:
        ObjectStoreReader.setProviders(null)
    }

    def 'both capabilities are optional and default to null'() {
        given: 'a provider that handles the scheme but implements neither operation'
        def provider = new ObjectStoreReader() { boolean canHandle(String s) { true } }

        expect: 'listWithMeta -> the caller falls back to the hierarchical NIO walk'
        provider.listWithMeta(Paths.get('/tmp/x')) == null
        and: 'readRange -> the caller falls back to a full read'
        provider.readRange(Paths.get('/tmp/x'), 0, 16384) == null
    }

    def 'discovery is not memoized, so a provider registered later is seen'() {
        given: 'no provider injected'
        ObjectStoreReader.setProviders(null)

        when: 'a look-up happens before the plugins are up and finds nothing'
        def first = ObjectStoreReader.getProviders()

        then: 'nothing is cached -- plugins start lazily per scheme, so the set still grows'
        first.isEmpty()
        cachedProviders() == null

        when: 'the provider becomes available'
        def provider = Stub(ObjectStoreReader) { canHandle('s3') >> true }
        ObjectStoreReader.setProviders([provider])

        then: 'it resolves; the earlier empty answer did not latch'
        ObjectStoreReader.getProviders() == [provider]

        cleanup:
        ObjectStoreReader.setProviders(null)
    }

    def 'readRange rejects a non-positive length, so no provider can issue an inverted range'() {
        given: 'a provider whose ranged read would be reached only if the check passed'
        def provider = new FakeReader()

        when: 'a valid window'
        def bytes = provider.readRange(Paths.get('/tmp/x'), 100, 3)
        then:
        bytes == [1,2,3] as byte[]

        when: 'len is zero -- `bytes=100-99` is inverted, and S3/Azure answer with the WHOLE object'
        provider.readRange(Paths.get('/tmp/x'), 100, 0)
        then:
        thrown(IllegalArgumentException)

        when: 'len is negative'
        provider.readRange(Paths.get('/tmp/x'), 100, -1)
        then:
        thrown(IllegalArgumentException)

        when: 'the offset is negative'
        provider.readRange(Paths.get('/tmp/x'), -1, 16384)
        then:
        thrown(IllegalArgumentException)
    }

    def 'normalizePrefix and relativize are shared so every cloud derives the same members'() {
        expect: 'a prefix ends with exactly one slash; the bucket root stays empty'
        ObjectStoreReader.normalizePrefix('dir') == 'dir/'
        ObjectStoreReader.normalizePrefix('dir/') == 'dir/'
        ObjectStoreReader.normalizePrefix('') == ''
        ObjectStoreReader.normalizePrefix(null) == ''

        and: 'members are relative to the prefix'
        ObjectStoreReader.relativize('dir/', 'dir/a.txt') == 'a.txt'
        ObjectStoreReader.relativize('dir/', 'dir/sub/b.txt') == 'sub/b.txt'
        ObjectStoreReader.relativize('', 'a.txt') == 'a.txt'

        and: 'a placeholder object and the prefix itself are skipped, not returned as members'
        ObjectStoreReader.relativize('dir/', 'dir/') == null
        ObjectStoreReader.relativize('dir/', 'dir/sub/') == null
    }

    /** Read the injection field directly: that discovery is not cached is not observable otherwise. */
    private static List<ObjectStoreReader> cachedProviders() {
        final field = ObjectStoreReader.getDeclaredField('providers')
        field.setAccessible(true)
        return (List<ObjectStoreReader>) field.get(null)
    }

}
