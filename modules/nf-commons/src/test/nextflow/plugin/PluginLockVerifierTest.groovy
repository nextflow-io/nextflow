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

package nextflow.plugin

import java.nio.file.Files
import java.nio.file.Path

import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PluginLockVerifierTest extends Specification {

    private Path zipWith(String content) {
        final zip = Files.createTempFile('plugin', '.zip')
        Files.write(zip, content.getBytes('UTF-8'))
        return zip
    }

    /** Write a lock file holding the given fqid->sha512 entries and return its path */
    private Path lockFileWith(Map<String,String> entries) {
        final path = Files.createTempFile('plugins', '.lock')
        final lock = new PluginLockFile()
        entries.each { k, v -> lock.addEntry(k, new PluginLockFile.Entry(v)) }
        lock.write(path)
        return path
    }

    @Unroll
    def 'should resolve lock mode from env [#VALUE]' () {
        given:
        SysEnv.push(VALUE != null ? [NXF_PLUGINS_LOCK_MODE: VALUE] : [:])

        expect:
        PluginLockVerifier.getMode() == EXPECTED

        cleanup:
        SysEnv.pop()

        where:
        VALUE       | EXPECTED
        null        | PluginLockVerifier.Mode.WARN
        'warn'      | PluginLockVerifier.Mode.WARN
        'WARN'      | PluginLockVerifier.Mode.WARN
        'strict'    | PluginLockVerifier.Mode.STRICT
        'STRICT'    | PluginLockVerifier.Mode.STRICT
        'off'       | PluginLockVerifier.Mode.OFF
        'OFF'       | PluginLockVerifier.Mode.OFF
        'nonsense'  | PluginLockVerifier.Mode.WARN
    }

    // ---------------------------------------------------------------------------
    // dormant / enabled

    def 'should be dormant when no lock file exists' () {
        given:
        SysEnv.push([NXF_PLUGINS_LOCK_MODE: 'strict'])
        and:
        def zip = zipWith('hello plugin')
        def verifier = new PluginLockVerifier(Path.of('/no/such/plugins.lock'))

        expect:
        !verifier.isEnabled()

        when:
        verifier.verify('nf-foo@1.0.0', zip)
        then:
        noExceptionThrown()

        cleanup:
        SysEnv.pop()
        Files.deleteIfExists(zip)
    }

    def 'should be enabled when an (even empty) lock file exists' () {
        given:
        def path = Files.createTempFile('plugins', '.lock')
        path.text = ''

        expect:
        new PluginLockVerifier(path).isEnabled()

        cleanup:
        Files.deleteIfExists(path)
    }

    // ---------------------------------------------------------------------------
    // trust-on-first-use pinning

    def 'should pin a coordinate missing from the lock on first sight (TOFU)' () {
        given:
        SysEnv.push([NXF_PLUGINS_LOCK_MODE: 'strict'])
        and:
        def zip = zipWith('hello plugin')
        final sha = PluginLockVerifier.sha512(zip)
        and:
        // an empty but present lock file -> feature enabled, nothing pinned yet
        def lockPath = lockFileWith([:])
        def verifier = new PluginLockVerifier(lockPath)

        when:
        verifier.verify('nf-foo@1.0.0', zip)

        then:
        noExceptionThrown()
        and:
        // the entry was appended to the lock file on disk
        def written = PluginLockFile.read(lockPath)
        written.getEntry('nf-foo@1.0.0').sha512 == sha

        cleanup:
        SysEnv.pop()
        Files.deleteIfExists(zip)
        Files.deleteIfExists(lockPath)
    }

    def 'should not pin when dormant' () {
        given:
        def zip = zipWith('hello plugin')
        def verifier = new PluginLockVerifier(Path.of('/no/such/plugins.lock'))

        when:
        verifier.verify('nf-foo@1.0.0', zip)

        then:
        !verifier.isEnabled()
        noExceptionThrown()

        cleanup:
        Files.deleteIfExists(zip)
    }

    // ---------------------------------------------------------------------------
    // verification of an already-locked coordinate

    def 'should pass verification when the artifact matches the lock' () {
        given:
        SysEnv.push([NXF_PLUGINS_LOCK_MODE: 'strict'])
        and:
        def zip = zipWith('hello plugin')
        final sha = PluginLockVerifier.sha512(zip)
        def lockPath = lockFileWith(['nf-foo@1.0.0': sha])
        def verifier = new PluginLockVerifier(lockPath)

        when:
        verifier.verify('nf-foo@1.0.0', zip)

        then:
        noExceptionThrown()

        cleanup:
        SysEnv.pop()
        Files.deleteIfExists(zip)
        Files.deleteIfExists(lockPath)
    }

    def 'strict mode should abort on a checksum mismatch' () {
        given:
        SysEnv.push([NXF_PLUGINS_LOCK_MODE: 'strict'])
        and:
        def zip = zipWith('hello plugin')
        def lockPath = lockFileWith(['nf-foo@1.0.0': 'deadbeef'])
        def verifier = new PluginLockVerifier(lockPath)

        when:
        verifier.verify('nf-foo@1.0.0', zip)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('checksum does not match')

        cleanup:
        SysEnv.pop()
        Files.deleteIfExists(zip)
        Files.deleteIfExists(lockPath)
    }

    def 'warn mode should not abort on a checksum mismatch' () {
        given:
        SysEnv.push([NXF_PLUGINS_LOCK_MODE: 'warn'])
        and:
        def zip = zipWith('hello plugin')
        def lockPath = lockFileWith(['nf-foo@1.0.0': 'deadbeef'])
        def verifier = new PluginLockVerifier(lockPath)

        when:
        verifier.verify('nf-foo@1.0.0', zip)
        // second call to exercise the "log once" path
        verifier.verify('nf-foo@1.0.0', zip)

        then:
        noExceptionThrown()

        cleanup:
        SysEnv.pop()
        Files.deleteIfExists(zip)
        Files.deleteIfExists(lockPath)
    }

    def 'off mode should ignore a checksum mismatch' () {
        given:
        SysEnv.push([NXF_PLUGINS_LOCK_MODE: 'off'])
        and:
        def zip = zipWith('hello plugin')
        def lockPath = lockFileWith(['nf-foo@1.0.0': 'deadbeef'])
        def verifier = new PluginLockVerifier(lockPath)

        when:
        verifier.verify('nf-foo@1.0.0', zip)

        then:
        noExceptionThrown()

        cleanup:
        SysEnv.pop()
        Files.deleteIfExists(zip)
        Files.deleteIfExists(lockPath)
    }

    // ---------------------------------------------------------------------------
    // archive-absent (e.g. a cache extracted before this feature): never abort

    def 'strict mode should NOT abort when a locked artifact is unavailable' () {
        given:
        SysEnv.push([NXF_PLUGINS_LOCK_MODE: 'strict'])
        and:
        final missing = Path.of('/no/such/nf-foo-1.0.0.zip')
        def lockPath = lockFileWith(['nf-foo@1.0.0': 'abc'])
        def verifier = new PluginLockVerifier(lockPath)

        when:
        verifier.verify('nf-foo@1.0.0', missing)

        then:
        noExceptionThrown()

        cleanup:
        SysEnv.pop()
        Files.deleteIfExists(lockPath)
    }

    def 'sha512 should compute a 128-char lowercase hex digest' () {
        given:
        def zip = zipWith('some bytes')

        when:
        final sha = PluginLockVerifier.sha512(zip)

        then:
        sha.length() == 128
        sha ==~ /[0-9a-f]{128}/

        cleanup:
        Files.deleteIfExists(zip)
    }
}
