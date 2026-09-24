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
        @Override protected boolean isEnabled(Map config) { false }
        @Override protected CacheDB newInstance(UUID uniqueId, String runName, Path home) { null }
    }

    /** A realistic adopter: enabled by its own config key, and nothing else. */
    static class ConfiguredFactory extends CacheFactory {
        @Override protected boolean isEnabled(Map config) { config?.mycache == true }
        @Override protected CacheDB newInstance(UUID uniqueId, String runName, Path home) { null }
    }

    def 'a factory that does not override isEnabled is selected, as before'() {
        given:
        def factory = new AlwaysFactory()

        expect:
        CacheFactory.select([factory], [:]) is factory
    }

    def 'a declining factory is skipped and the next one serves'() {
        given:
        def declining = new DecliningFactory()
        def serving = new AlwaysFactory()

        expect: 'the declining factory is first in priority order, and still not chosen'
        CacheFactory.select([declining, serving], [:]) is serving
    }

    def 'priority order decides among the factories that claim the session'() {
        given:
        def first = new AlwaysFactory()
        def second = new AlwaysFactory()

        expect:
        CacheFactory.select([first, second], [:]) is first
    }

    def 'aborts when every registered factory declines, naming them'() {
        when:
        CacheFactory.select([new DecliningFactory(), new DecliningFactory()], [:])

        then:
        def e = thrown(IllegalStateException)
        e.message.contains('Unable to find an enabled Nextflow cache factory')
        e.message.contains(DecliningFactory.getName())
    }

    def 'the config decides, and only the config'() {
        given:
        def configured = new ConfiguredFactory()
        def fallback = new AlwaysFactory()

        expect: 'its key is set -- it serves'
        CacheFactory.select([configured, fallback], [mycache: true]) is configured

        and: 'its key is absent -- the next one serves'
        CacheFactory.select([configured, fallback], [something: 'else']) is fallback

        and: 'outside a session there is no config at all, which is not a claim'
        CacheFactory.select([configured, fallback], null) is fallback
    }

    def 'the same factory is chosen at init and at cleanup'() {
        given: '''Session.init creates the cache and Session.cleanup creates it AGAIN, after the
              session was destroyed. Between the two, newInstance is allowed to write into the
              session -- workDir, resumeMode. Reading only the config is what makes the second
              selection reach the same factory as the first, so the records are read back by the
              backend that wrote them'''
        def config = [mycache: true]
        def factories = [new ConfiguredFactory(), new AlwaysFactory()]

        when:
        def atInit = CacheFactory.select(factories, config)
        and: 'whatever the session went through in between, the config is the same object'
        def atCleanup = CacheFactory.select(factories, config)

        then:
        atCleanup.is(atInit)
    }

    def 'aborts when no factory is registered at all'() {
        when:
        CacheFactory.select(EMPTY, [:])

        then:
        def e = thrown(IllegalStateException)
        e.message == 'Unable to find Nextflow cache factory'

        where:
        EMPTY << [ [] as List<CacheFactory>, null ]
    }

}
