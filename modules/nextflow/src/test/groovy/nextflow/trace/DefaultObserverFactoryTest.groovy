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

package nextflow.trace

import nextflow.Session
import nextflow.cache.CacheDB
import nextflow.cache.CacheStore
import nextflow.exception.AbortOperationException
import spock.lang.Specification

class DefaultObserverFactoryTest extends Specification {

    def 'should create cleanup observer only with eager cleanup and the default cache' () {
        given:
        def session = Mock(Session) {
            getConfig() >> [cleanup: CLEANUP]
            getCache() >> CACHE
        }

        when:
        def observers = new DefaultObserverFactory().create(session)
        then:
        observers.any { o -> o instanceof CleanupObserver } == EXPECTED

        where:
        CLEANUP     | CACHE                                 | EXPECTED
        null        | null                                  | false
        true        | null                                  | false
        'eager'     | null                                  | true
        'eager'     | new CacheDB(Mock(CacheStore))         | true
        'eager'     | new CustomCacheDB(Mock(CacheStore))   | false
    }

    def 'should reject an invalid cleanup option' () {
        given:
        def session = Mock(Session) {
            getConfig() >> [cleanup: 'lazy']
        }

        when:
        new DefaultObserverFactory().create(session)
        then:
        thrown(AbortOperationException)
    }

}

class CustomCacheDB extends CacheDB {

    CustomCacheDB(CacheStore store) {
        super(store)
    }

}
