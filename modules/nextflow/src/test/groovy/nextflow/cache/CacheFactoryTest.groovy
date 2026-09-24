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

import java.nio.file.Path

import spock.lang.Specification

/**
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class CacheFactoryTest extends Specification {

    /** A factory that keeps the inherited isEnabled(), i.e. an existing one. */
    static class AlwaysFactory extends CacheFactory {
        @Override protected CacheDB newInstance(UUID uniqueId, String runName, Path home) { null }
    }

    /** A factory that declines, as one reading its own configuration would. */
    static class DecliningFactory extends CacheFactory {
        @Override protected boolean isEnabled() { false }
        @Override protected CacheDB newInstance(UUID uniqueId, String runName, Path home) { null }
    }

    def 'a factory that does not override isEnabled is selected, as before'() {
        given:
        def factory = new AlwaysFactory()

        expect:
        CacheFactory.select([factory]) is factory
    }

    def 'a declining factory is skipped and the next one serves'() {
        given:
        def declining = new DecliningFactory()
        def serving = new AlwaysFactory()

        expect: 'the declining factory is first in priority order, and still not chosen'
        CacheFactory.select([declining, serving]) is serving
    }

    def 'priority order decides among the factories that claim the session'() {
        given:
        def first = new AlwaysFactory()
        def second = new AlwaysFactory()

        expect:
        CacheFactory.select([first, second]) is first
    }

    def 'aborts when every registered factory declines, naming them'() {
        when:
        CacheFactory.select([new DecliningFactory(), new DecliningFactory()])

        then:
        def e = thrown(IllegalStateException)
        e.message.contains('Unable to find an enabled Nextflow cache factory')
        e.message.contains(DecliningFactory.getName())
    }

    def 'aborts when no factory is registered at all'() {
        when:
        CacheFactory.select(EMPTY)

        then:
        def e = thrown(IllegalStateException)
        e.message == 'Unable to find Nextflow cache factory'

        where:
        EMPTY << [ [] as List<CacheFactory>, null ]
    }

}
