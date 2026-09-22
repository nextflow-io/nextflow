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

import java.nio.file.FileSystem
import java.nio.file.Path
import java.nio.file.spi.FileSystemProvider

import spock.lang.Specification

class AtomicLockProviderTest extends Specification {

    def cleanup() {
        AtomicLockProvider.setProviders(null)
    }

    private Path schemePath(String scheme) {
        Stub(Path) {
            getFileSystem() >> Stub(FileSystem) {
                provider() >> Stub(FileSystemProvider) { getScheme() >> scheme }
            }
        }
    }

    def 'lookup resolves the provider by scheme'() {
        given:
        def provider = Stub(AtomicLockProvider) { canHandle('s3') >> true }
        AtomicLockProvider.setProviders([provider])

        expect:
        AtomicLockProvider.lookup(schemePath('s3')).is(provider)
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
        def provider = Stub(AtomicLockProvider) { canHandle('s3') >> true }
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
        AtomicLockProvider.lookup(schemePath('zz'))
        then:
        thrown(IllegalStateException)
    }

}
