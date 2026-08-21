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

    /** Create an extracted-plugin-like directory tree with the given relative-path->content files */
    private Path pluginDir(Map<String,String> files) {
        final dir = Files.createTempDirectory('plugin')
        files.each { rel, content ->
            final f = dir.resolve(rel)
            Files.createDirectories(f.parent)
            Files.write(f, content.getBytes('UTF-8'))
        }
        return dir
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
    // tree hashing

    def 'sha512Tree should be deterministic and content-sensitive' () {
        given:
        def a = pluginDir(['classes/A.class': 'aaa', 'META-INF/MANIFEST.MF': 'mmm'])
        def b = pluginDir(['classes/A.class': 'aaa', 'META-INF/MANIFEST.MF': 'mmm'])
        def c = pluginDir(['classes/A.class': 'aaa', 'META-INF/MANIFEST.MF': 'CHANGED'])

        expect:
        // identical content in different dirs -> identical hash
        PluginLockVerifier.sha512Tree(a) == PluginLockVerifier.sha512Tree(b)
        // any content change -> different hash
        PluginLockVerifier.sha512Tree(a) != PluginLockVerifier.sha512Tree(c)
        and:
        PluginLockVerifier.sha512Tree(a).length() == 128
        PluginLockVerifier.sha512Tree(a) ==~ /[0-9a-f]{128}/

        cleanup:
        [a, b, c].each { it.deleteDir() }
    }

    def 'sha512Tree should skip non-regular files' () {
        given:
        def dir = pluginDir(['classes/A.class': 'aaa'])
        final expected = PluginLockVerifier.sha512Tree(dir)
        and:
        // a symlink to a directory cannot be read as a byte stream: it must be skipped
        Files.createSymbolicLink(dir.resolve('link'), dir.resolve('classes'))

        expect:
        PluginLockVerifier.sha512Tree(dir) == expected

        cleanup:
        dir.deleteDir()
    }

    // ---------------------------------------------------------------------------
    // dormant / enabled

    def 'should be dormant when no lock file exists' () {
        given:
        SysEnv.push([NXF_PLUGINS_LOCK_MODE: 'strict'])
        and:
        def dir = pluginDir(['classes/A.class': 'aaa'])
        def verifier = new PluginLockVerifier(Path.of('/no/such/plugins.lock'))

        expect:
        !verifier.isEnabled()
        and:
        // no lock file location at all ie. the plugin system was never given a project dir
        !new PluginLockVerifier(null).isEnabled()

        when:
        verifier.verify('nf-foo@1.0.0', dir)
        then:
        noExceptionThrown()

        cleanup:
        SysEnv.pop()
        dir.deleteDir()
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
        def dir = pluginDir(['classes/A.class': 'hello'])
        final sha = PluginLockVerifier.sha512Tree(dir)
        and:
        // an empty but present lock file -> feature enabled, nothing pinned yet
        def lockPath = lockFileWith([:])
        def verifier = new PluginLockVerifier(lockPath)

        when:
        verifier.verify('nf-foo@1.0.0', dir)

        then:
        noExceptionThrown()
        and:
        // the entry was appended to the lock file on disk
        def written = PluginLockFile.read(lockPath)
        written.getEntry('nf-foo@1.0.0').sha512 == sha

        cleanup:
        SysEnv.pop()
        dir.deleteDir()
        Files.deleteIfExists(lockPath)
    }

    // ---------------------------------------------------------------------------
    // verification of an already-locked coordinate

    def 'should pass verification when the plugin matches the lock' () {
        given:
        SysEnv.push([NXF_PLUGINS_LOCK_MODE: 'strict'])
        and:
        def dir = pluginDir(['classes/A.class': 'hello'])
        final sha = PluginLockVerifier.sha512Tree(dir)
        def lockPath = lockFileWith(['nf-foo@1.0.0': sha])
        def verifier = new PluginLockVerifier(lockPath)

        when:
        verifier.verify('nf-foo@1.0.0', dir)

        then:
        noExceptionThrown()

        cleanup:
        SysEnv.pop()
        dir.deleteDir()
        Files.deleteIfExists(lockPath)
    }

    def 'strict mode should abort when the extracted plugin does not match' () {
        given:
        SysEnv.push([NXF_PLUGINS_LOCK_MODE: 'strict'])
        and:
        def dir = pluginDir(['classes/A.class': 'poisoned'])
        def lockPath = lockFileWith(['nf-foo@1.0.0': 'deadbeef'])
        def verifier = new PluginLockVerifier(lockPath)

        when:
        verifier.verify('nf-foo@1.0.0', dir)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('does not match')

        cleanup:
        SysEnv.pop()
        dir.deleteDir()
        Files.deleteIfExists(lockPath)
    }

    def 'warn mode should not abort on a mismatch' () {
        given:
        SysEnv.push([NXF_PLUGINS_LOCK_MODE: 'warn'])
        and:
        def dir = pluginDir(['classes/A.class': 'poisoned'])
        def lockPath = lockFileWith(['nf-foo@1.0.0': 'deadbeef'])
        def verifier = new PluginLockVerifier(lockPath)

        when:
        verifier.verify('nf-foo@1.0.0', dir)
        // second call to exercise the "log once" path
        verifier.verify('nf-foo@1.0.0', dir)

        then:
        noExceptionThrown()

        cleanup:
        SysEnv.pop()
        dir.deleteDir()
        Files.deleteIfExists(lockPath)
    }

    def 'off mode should disable the verifier altogether' () {
        given:
        SysEnv.push([NXF_PLUGINS_LOCK_MODE: 'off'])
        and:
        def dir = pluginDir(['classes/A.class': 'poisoned'])
        def lockPath = lockFileWith(['nf-foo@1.0.0': 'deadbeef'])
        def verifier = new PluginLockVerifier(lockPath)
        final before = lockPath.text

        expect:
        !verifier.isEnabled()

        when:
        // both a mismatch and an unlocked coordinate are no-ops
        verifier.verify('nf-foo@1.0.0', dir)
        verifier.verify('nf-bar@1.0.0', dir)

        then:
        noExceptionThrown()
        and:
        // nothing was pinned ie. the user working tree is left untouched
        lockPath.text == before

        cleanup:
        SysEnv.pop()
        dir.deleteDir()
        Files.deleteIfExists(lockPath)
    }

    def 'off mode should not even read a malformed lock file' () {
        given:
        SysEnv.push([NXF_PLUGINS_LOCK_MODE: 'off'])
        and:
        def lockPath = Files.createTempFile('plugins', '.lock')
        lockPath.text = '{ this is not valid json ]'

        when:
        def verifier = new PluginLockVerifier(lockPath)

        then:
        noExceptionThrown()
        !verifier.isEnabled()

        cleanup:
        SysEnv.pop()
        Files.deleteIfExists(lockPath)
    }

    def 'should warn instead of aborting when the lock file cannot be written' () {
        given:
        SysEnv.push([NXF_PLUGINS_LOCK_MODE: 'strict'])
        and:
        def dir = pluginDir(['classes/A.class': 'hello'])
        def folder = Files.createTempDirectory('test')
        def lockPath = folder.resolve('plugins.lock')
        lockPath.text = ''
        def verifier = new PluginLockVerifier(lockPath)
        and:
        // make the lock file location unwritable after the verifier has been created
        folder.deleteDir()

        when:
        verifier.verify('nf-foo@1.0.0', dir)

        then:
        noExceptionThrown()

        cleanup:
        SysEnv.pop()
        dir.deleteDir()
        folder?.deleteDir()
    }

    def 'should catch cache poisoning regardless of directory ownership' () {
        given:
        // pin the good tree, then tamper an extracted file in place (the shared-cache attack)
        SysEnv.push([NXF_PLUGINS_LOCK_MODE: 'strict'])
        and:
        def dir = pluginDir(['classes/A.class': 'good', 'lib/x.jar': 'jar'])
        final good = PluginLockVerifier.sha512Tree(dir)
        def lockPath = lockFileWith(['nf-foo@1.0.0': good])
        def verifier = new PluginLockVerifier(lockPath)
        and:
        // attacker edits an already-extracted file
        Files.write(dir.resolve('classes/A.class'), 'evil'.getBytes('UTF-8'))

        when:
        verifier.verify('nf-foo@1.0.0', dir)

        then:
        thrown(AbortOperationException)

        cleanup:
        SysEnv.pop()
        dir.deleteDir()
        Files.deleteIfExists(lockPath)
    }
}
