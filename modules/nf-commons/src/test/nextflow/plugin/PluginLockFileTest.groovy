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

import nextflow.exception.AbortOperationException
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PluginLockFileTest extends Specification {

    def 'should round-trip write and read' () {
        given:
        def path = Files.createTempFile('plugins', '.lock')
        and:
        final lock = new PluginLockFile()
        lock.version = 1
        lock.addEntry('nf-amazon@2.0.0', new PluginLockFile.Entry('sha512-aaa'))
        lock.addEntry('nf-google@1.5.0', new PluginLockFile.Entry('sha512-ggg'))

        when:
        lock.write(path)
        and:
        final copy = PluginLockFile.read(path)

        then:
        copy.version == 1
        copy.getEntry('nf-amazon@2.0.0') == new PluginLockFile.Entry('sha512-aaa')
        copy.getEntry('nf-google@1.5.0') == new PluginLockFile.Entry('sha512-ggg')
        copy.entries.size() == 2
        and:
        copy == lock

        cleanup:
        Files.deleteIfExists(path)
    }

    def 'should write pretty json with stable key order' () {
        given:
        def path = Files.createTempFile('plugins', '.lock')
        and:
        final lock = new PluginLockFile()
        lock.version = 1
        // add out of order
        lock.addEntry('nf-zeta@1.0.0', new PluginLockFile.Entry('sha512-z'))
        lock.addEntry('nf-alpha@1.0.0', new PluginLockFile.Entry('sha512-a'))

        when:
        lock.write(path)
        final text = path.text

        then:
        // pretty printed (contains newlines and indentation)
        text.contains('\n')
        // keys are sorted: alpha appears before zeta
        text.indexOf('nf-alpha@1.0.0') < text.indexOf('nf-zeta@1.0.0')

        cleanup:
        Files.deleteIfExists(path)
    }

    def 'should return empty for missing file' () {
        given:
        final path = Path.of('/no/such/plugins.lock')

        when:
        final lock = PluginLockFile.read(path)

        then:
        lock != null
        lock.isEmpty()
        lock.entries.isEmpty()
    }

    def 'should return empty for a blank (touched) file' () {
        given:
        def path = Files.createTempFile('plugins', '.lock')
        path.text = '   \n'

        when:
        final lock = PluginLockFile.read(path)

        then:
        lock != null
        lock.isEmpty()

        cleanup:
        Files.deleteIfExists(path)
    }

    def 'should throw for malformed file' () {
        given:
        def path = Files.createTempFile('plugins', '.lock')
        path.text = '{ this is not valid json ]'

        when:
        PluginLockFile.read(path)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('malformed JSON')

        cleanup:
        Files.deleteIfExists(path)
    }

    def 'should throw for a lock file created by a newer version' () {
        given:
        def path = Files.createTempFile('plugins', '.lock')
        path.text = '{ "version": 2, "plugins": {} }'

        when:
        PluginLockFile.read(path)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('requires a newer version of Nextflow')

        cleanup:
        Files.deleteIfExists(path)
    }

    def 'should replace the target file leaving no temporary file behind' () {
        given:
        def folder = Files.createTempDirectory('test')
        def path = folder.resolve('plugins.lock')
        path.text = '{ "version": 1, "plugins": { "nf-old@1.0.0": { "sha512": "sha512-old" } } }'

        when:
        final lock = new PluginLockFile()
        lock.addEntry('nf-new@1.0.0', new PluginLockFile.Entry('sha512-new'))
        lock.write(path)

        then:
        // the previous content is fully replaced and no leftover temp file remains
        PluginLockFile.read(path).entries.keySet() == ['nf-new@1.0.0'] as Set
        and:
        Files.list(folder).collect { it.fileName.toString() } == ['plugins.lock']

        cleanup:
        folder?.deleteDir()
    }

    def 'should lookup an entry by id and version' () {
        given:
        final lock = new PluginLockFile()
        lock.addEntry('nf-amazon@2.0.0', new PluginLockFile.Entry('sha512-aaa'))

        expect:
        lock.getEntry('nf-amazon@2.0.0') == new PluginLockFile.Entry('sha512-aaa')
        lock.getEntry('nf-unknown@9.9.9') == null
    }

    def 'should add and update an entry' () {
        given:
        final lock = new PluginLockFile()

        when:
        lock.addEntry('nf-amazon@2.0.0', new PluginLockFile.Entry('sha512-old'))
        then:
        lock.getEntry('nf-amazon@2.0.0').sha512 == 'sha512-old'
        lock.entries.size() == 1

        when:
        lock.addEntry('nf-amazon@2.0.0', new PluginLockFile.Entry('sha512-new'))
        then:
        lock.getEntry('nf-amazon@2.0.0').sha512 == 'sha512-new'
        lock.entries.size() == 1
    }

}
