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

class CloudCacheStoreTest extends Specification {

    def 'should get and put hash index entries' () {
        given:
        def base = Files.createTempDirectory('cloudcache')
        // note: open() is not needed - the hash-index methods are independent of
        // the run-index writer and create their own parent directory on write
        def store = new CloudCacheStore(UUID.randomUUID(), 'run_1', base)
        and:
        def content = CacheHelper.hasher('CONTENT').hash()
        def finalHash = CacheHelper.hasher('FINAL').hash()

        expect:
        store.getSuccessfulHash(content) == null

        when:
        store.putSuccessfulHash(content, finalHash)
        then:
        store.getSuccessfulHash(content) == finalHash

        when:
        store.deleteSuccessfulHash(content)
        then:
        store.getSuccessfulHash(content) == null

        cleanup:
        store?.close()
        base?.deleteDir()
    }

    def 'hash index is scoped per session (uniqueId)' () {
        given:
        def base = Files.createTempDirectory('cloudcache')
        def content = CacheHelper.hasher('CONTENT').hash()
        def finalHash = CacheHelper.hasher('FINAL').hash()
        and: 'session A writes a pointer'
        def storeA = new CloudCacheStore(UUID.randomUUID(), 'run_A', base)
        storeA.putSuccessfulHash(content, finalHash)
        and: 'a different session B over the same base'
        def storeB = new CloudCacheStore(UUID.randomUUID(), 'run_B', base)

        expect: 'session B does not see session A pointer'
        storeB.getSuccessfulHash(content) == null

        cleanup:
        storeA?.close(); storeB?.close()
        base?.deleteDir()
    }
}
