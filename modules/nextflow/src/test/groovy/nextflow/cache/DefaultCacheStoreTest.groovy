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

package nextflow.cache

import java.nio.file.Files

import nextflow.util.CacheHelper
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DefaultCacheStoreTest extends Specification {

    def 'should get and put cache entries' () {
        given:
        def folder = Files.createTempDirectory('test')
        def uuid = UUID.randomUUID()
        def runName = 'test_1'
        and:
        def store = new DefaultCacheStore(uuid, runName, folder); store.open()

        and:
        def key1 = CacheHelper.hasher('ONE').hash()
        def key2 = CacheHelper.hasher('TWO').hash()
        def value = "Hello world"

        when:
        store.putEntry(key1, value.bytes)
        then:
        new String(store.getEntry(key1)) == value
        and:
        store.getEntry(key2) == null

        when:
        store.deleteEntry(key1)
        then:
        store.getEntry(key1) == null

        cleanup:
        store?.close()
        folder?.deleteDir()
    }

    def 'should get and put hash index entries' () {
        given:
        def folder = Files.createTempDirectory('test')
        def uuid = UUID.randomUUID()
        def store = new DefaultCacheStore(uuid, 'test_1', folder); store.open()
        and:
        def content = CacheHelper.hasher('CONTENT').hash()
        def finalHash = CacheHelper.hasher('FINAL').hash()

        expect: 'absent key returns null'
        store.getHashIndex(content) == null

        when:
        store.putHashIndex(content, finalHash)
        then: 'round-trips'
        store.getHashIndex(content) == finalHash

        and: 'index key does not collide with an entry key of the same hash'
        store.getEntry(content) == null
        store.putEntry(content, 'ENTRY'.bytes)
        new String(store.getEntry(content)) == 'ENTRY'
        store.getHashIndex(content) == finalHash

        when: 'the cache is dropped'
        store.close()
        store.drop()
        and: 'reopened on the same session location'
        store = new DefaultCacheStore(uuid, 'test_1', folder); store.open()
        then: 'the index entry is gone'
        store.getHashIndex(content) == null

        cleanup:
        store?.close()
        folder?.deleteDir()
    }

}
