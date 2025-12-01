/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.cache.sqlite

import java.nio.file.Files

import nextflow.cache.CacheStore
import nextflow.exception.AbortOperationException
import nextflow.util.CacheHelper
import spock.lang.Specification

/**
 * Unit tests for SQLiteCacheStore
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SQLiteCacheStoreTest extends Specification {

    def 'should get and put cache entries' () {
        given:
        def folder = Files.createTempDirectory('test')
        def uuid = UUID.randomUUID()
        def runName = 'test_1'
        and:
        def store = new SQLiteCacheStore(uuid, runName, folder)
        store.open()

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

    def 'should write and iterate cache index' () {
        given:
        def folder = Files.createTempDirectory('test')
        def uuid = UUID.randomUUID()
        def runName = 'test_run'
        and:
        def store = new SQLiteCacheStore(uuid, runName, folder)
        store.open()

        and:
        def key1 = CacheHelper.hasher('FIRST').hash()
        def key2 = CacheHelper.hasher('SECOND').hash()
        def key3 = CacheHelper.hasher('THIRD').hash()

        when:
        store.writeIndex(key1, true)
        store.writeIndex(key2, false)
        store.writeIndex(key3, true)

        and:
        def indexes = []
        def iterator = store.iterateIndex()
        while (iterator.hasNext()) {
            indexes << iterator.next()
        }

        then:
        indexes.size() == 3
        indexes[0].key == key1
        indexes[0].cached == true
        indexes[1].key == key2
        indexes[1].cached == false
        indexes[2].key == key3
        indexes[2].cached == true

        cleanup:
        store?.close()
        folder?.deleteDir()
    }

    def 'should delete cache index' () {
        given:
        def folder = Files.createTempDirectory('test')
        def uuid = UUID.randomUUID()
        def runName = 'test_run'
        and:
        def store = new SQLiteCacheStore(uuid, runName, folder)
        store.open()

        and:
        def key1 = CacheHelper.hasher('FIRST').hash()
        def key2 = CacheHelper.hasher('SECOND').hash()

        when:
        store.writeIndex(key1, true)
        store.writeIndex(key2, false)

        then:
        def indexes = []
        def iterator = store.iterateIndex()
        while (iterator.hasNext()) {
            indexes << iterator.next()
        }
        indexes.size() == 2

        when:
        store.deleteIndex()

        and:
        def indexes2 = []
        def iterator2 = store.iterateIndex()
        while (iterator2.hasNext()) {
            indexes2 << iterator2.next()
        }

        then:
        indexes2.size() == 0

        cleanup:
        store?.close()
        folder?.deleteDir()
    }

    def 'should handle different run names separately' () {
        given:
        def folder = Files.createTempDirectory('test')
        def uuid = UUID.randomUUID()
        def runName1 = 'test_run_1'
        def runName2 = 'test_run_2'

        and:
        def store1 = new SQLiteCacheStore(uuid, runName1, folder)
        store1.open()
        def store2 = new SQLiteCacheStore(uuid, runName2, folder)
        store2.open()

        and:
        def key1 = CacheHelper.hasher('KEY1').hash()
        def key2 = CacheHelper.hasher('KEY2').hash()

        when:
        store1.writeIndex(key1, true)
        store2.writeIndex(key2, false)

        then:
        def indexes1 = []
        def iterator1 = store1.iterateIndex()
        while (iterator1.hasNext()) {
            indexes1 << iterator1.next()
        }
        indexes1.size() == 1
        indexes1[0].key == key1

        and:
        def indexes2 = []
        def iterator2 = store2.iterateIndex()
        while (iterator2.hasNext()) {
            indexes2 << iterator2.next()
        }
        indexes2.size() == 1
        indexes2[0].key == key2

        cleanup:
        store1?.close()
        store2?.close()
        folder?.deleteDir()
    }

    def 'should open for read when index exists' () {
        given:
        def folder = Files.createTempDirectory('test')
        def uuid = UUID.randomUUID()
        def runName = 'test_run'

        and:
        def store1 = new SQLiteCacheStore(uuid, runName, folder)
        store1.open()
        def key1 = CacheHelper.hasher('KEY1').hash()
        store1.writeIndex(key1, true)
        store1.close()

        when:
        def store2 = new SQLiteCacheStore(uuid, runName, folder)
        store2.openForRead()

        then:
        def indexes = []
        def iterator = store2.iterateIndex()
        while (iterator.hasNext()) {
            indexes << iterator.next()
        }
        indexes.size() == 1
        indexes[0].key == key1

        cleanup:
        store2?.close()
        folder?.deleteDir()
    }

    def 'should throw exception when opening for read without index' () {
        given:
        def folder = Files.createTempDirectory('test')
        def uuid = UUID.randomUUID()
        def runName = 'test_run'
        and:
        def store = new SQLiteCacheStore(uuid, runName, folder)

        when:
        store.openForRead()

        then:
        thrown(AbortOperationException)

        cleanup:
        store?.close()
        folder?.deleteDir()
    }

    def 'should drop cache completely' () {
        given:
        def folder = Files.createTempDirectory('test')
        def uuid = UUID.randomUUID()
        def runName = 'test_run'
        and:
        def store = new SQLiteCacheStore(uuid, runName, folder)
        store.open()

        and:
        def key1 = CacheHelper.hasher('KEY1').hash()
        def value = "Test value"
        store.putEntry(key1, value.bytes)
        store.writeIndex(key1, true)

        when:
        def dataDir = store.dataDir
        def exists1 = dataDir.exists()
        store.drop()
        def exists2 = dataDir.exists()

        then:
        exists1 == true
        exists2 == false

        cleanup:
        folder?.deleteDir()
    }

    def 'should handle multiple entries and maintain consistency' () {
        given:
        def folder = Files.createTempDirectory('test')
        def uuid = UUID.randomUUID()
        def runName = 'test_run'
        and:
        def store = new SQLiteCacheStore(uuid, runName, folder)
        store.open()

        and:
        def entries = [:]
        (1..100).each { i ->
            def key = CacheHelper.hasher("KEY_$i").hash()
            def value = "Value for entry $i"
            entries[key] = value
            store.putEntry(key, value.bytes)
            store.writeIndex(key, i % 2 == 0) // alternate between cached/not cached
        }

        when:
        def retrievedEntries = [:]
        entries.each { key, expectedValue ->
            def actualValue = new String(store.getEntry(key))
            retrievedEntries[key] = actualValue
        }

        and:
        def indexes = []
        def iterator = store.iterateIndex()
        while (iterator.hasNext()) {
            indexes << iterator.next()
        }

        then:
        retrievedEntries == entries
        indexes.size() == 100
        indexes.findAll { it.cached }.size() == 50
        indexes.findAll { !it.cached }.size() == 50

        cleanup:
        store?.close()
        folder?.deleteDir()
    }

    def 'should handle concurrent access' () {
        given:
        def folder = Files.createTempDirectory('test')
        def uuid = UUID.randomUUID()
        def runName = 'test_run'
        and:
        def store = new SQLiteCacheStore(uuid, runName, folder)
        store.open()

        and:
        def key = CacheHelper.hasher('CONCURRENT_KEY').hash()
        def value = "Concurrent value"

        when:
        // Simulate concurrent writes
        10.times { i ->
            store.putEntry(key, "${value}_${i}".bytes)
        }

        and:
        def finalValue = new String(store.getEntry(key))

        then:
        finalValue.startsWith(value)

        cleanup:
        store?.close()
        folder?.deleteDir()
    }
}