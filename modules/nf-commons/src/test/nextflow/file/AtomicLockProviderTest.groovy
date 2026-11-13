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

import spock.lang.Specification

class AtomicLockProviderTest extends Specification {

    def cleanup() {
        AtomicLockProvider.setProviders(null)
    }

    static class FakeLockProvider extends AtomicLockProvider {
        String scheme = 'file'
        boolean canHandle(String scheme) { scheme == this.scheme }
        boolean tryCreate(Path lockPath) { false }
    }

    def 'lookup resolves the provider by scheme'() {
        given:
        def provider = new FakeLockProvider()
        AtomicLockProvider.setProviders([new FakeLockProvider(scheme: 's3'), provider])

        expect:
        AtomicLockProvider.lookup(Path.of('/tmp/x')).is(provider)   // 'file' scheme
    }

    def 'discovery is not memoized, so a provider registered later is seen'() {
        given: 'no provider has been injected and the plugin system is not started'
        AtomicLockProvider.setProviders(null)

        when: 'the very first look-up happens too early and finds nothing'
        def first = AtomicLockProvider.getProviders()

        then: 'nothing is cached -- latching it would make lookup throw for a plugin that IS loaded'
        first.isEmpty()
        cachedProviders() == null

        when: 'the provider becomes available'
        def provider = new FakeLockProvider()
        AtomicLockProvider.setProviders([provider])

        then: 'it resolves; the earlier empty answer did not latch'
        AtomicLockProvider.getProviders() == [provider]
    }

    /** Read the injection field directly: that discovery is not cached is not observable otherwise. */
    private static List<AtomicLockProvider> cachedProviders() {
        final field = AtomicLockProvider.getDeclaredField('providers')
        field.setAccessible(true)
        return (List<AtomicLockProvider>) field.get(null)
    }

    def 'lookup throws when no provider serves the scheme'() {
        given:
        AtomicLockProvider.setProviders([])

        when:
        AtomicLockProvider.lookup(Path.of('/tmp/x'))
        then:
        thrown(IllegalStateException)
    }

}
