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

package nextflow.script.bundle

import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.Path

import nextflow.file.FileHelper
import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.lang.TempDir

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ResourcesBundleTest extends Specification {

    def LAST_MODIFIED = 1_000_000_000_000

    @TempDir
    Path folder

    def 'should scan bundle files' () {
        given:
        def bundlePath = folder.resolve('mod1'); bundlePath.mkdir()
        bundlePath.resolve('main.nf').text = "I'm the main file"
        bundlePath.resolve('this/that').mkdirs()
        Files.createFile(bundlePath.resolve('this/hola.txt'))
        Files.createFile(bundlePath.resolve('this/hello.txt'))
        Files.createFile(bundlePath.resolve('this/that/ciao.txt'))

        bundlePath.resolve('main.nf').setPermissions(6,4,4)
        bundlePath.resolve('this').setPermissions(7,5,5)
        bundlePath.resolve('this/that').setPermissions(7,5,5)
        bundlePath.resolve('this/hola.txt').setPermissions(6,4,4)
        bundlePath.resolve('this/hello.txt').setPermissions(6,4,4)
        bundlePath.resolve('this/that/ciao.txt').setPermissions(6,4,4)

        and:
        FileHelper.visitFiles([type:'any'], bundlePath, '**', { Path it -> it.setLastModified(LAST_MODIFIED) })

        when:
        def bundle = ResourcesBundle.scan(bundlePath)
        then:
        bundle
        bundle.hasEntries()
        and:
        bundle.getPaths() == [
                bundlePath.resolve('main.nf'),
                bundlePath.resolve('this'),
                bundlePath.resolve('this/hola.txt'),
                bundlePath.resolve('this/hello.txt'),
                bundlePath.resolve('this/that'),
                bundlePath.resolve('this/that/ciao.txt') ] as Set
        and:
        bundle.getEntries() == [
                'main.nf',
                'this',
                'this/hola.txt',
                'this/hello.txt',
                'this/that',
                'this/that/ciao.txt' ] as Set

        and:
        bundle.fingerprint() == '7404949a73f1707e39f641d9a54261ae'

    }

    def 'should get dockerfile' () {
        given:
        def dockerPath = folder.resolve('Dockerfile'); dockerPath.text = "I'm the main file"
        def bundlePath = folder.resolve('bundle')
        and:
        dockerPath.setLastModified(LAST_MODIFIED)
        dockerPath.setPermissions(6,4,4)
        when:
        def bundle = ResourcesBundle.scan(bundlePath)
        then:
        bundle.getDockerfile() == dockerPath
        and:
        bundle
        !bundle.hasEntries()
        and:
        bundle.fingerprint() == 'c5e75f95b68f5119debd926961ef6abd'

        when:
        // changing file permissions, change the fingerprint
        dockerPath.setPermissions(6,0,0)
        then:
        bundle.fingerprint() == '7b2200ff24230f76cea22e5eb15b1701'

    }

    def 'should get singularityfile' () {
        given:
        def singularPath = folder.resolve('Singularityfile'); singularPath.text = "I'm the main file"
        def bundlePath = folder.resolve('bundle')
        and:
        singularPath.setLastModified(LAST_MODIFIED)
        singularPath.setPermissions(6,4,4)
        when:
        def bundle = ResourcesBundle.scan(bundlePath)
        then:
        bundle.getSingularityfile() == singularPath
        and:
        bundle
        !bundle.hasEntries()
        and:
        bundle.fingerprint() == '6933e9238f3363c8e013a35715fa0540'

        when:
        // changing file permissions, change the fingerprint
        singularPath.setPermissions(6,0,0)
        then:
        bundle.fingerprint() == '3ffe7f16cd5ae17e6ba7485e01972b20'

    }

    def 'should check max file size'() {
        given:
        def root = folder.resolve('mod1'); root.mkdir()
        root.resolve('main.nf').text = "I'm the main file"

        when:
        ResourcesBundle.scan(root, [maxFileSize: MemoryUnit.of(5)])
        then:
        thrown(IllegalArgumentException)

    }

    def 'should check max bundle size'() {
        given:
        def root = folder.resolve('mod1'); root.mkdir()
        root.resolve('main.nf').text = "I'm the main file"

        when:
        ResourcesBundle.scan(root, [maxBundleSize: MemoryUnit.of(5)])
        then:
        def e = thrown(IllegalArgumentException)
        e.message == 'Module total size cannot exceed 5 B'

    }

    def 'should symlink not allowed'() {
        given:
        def root = folder.resolve('mod1'); root.mkdir()
        def main = root.resolve('main.nf'); main.text = "I'm the main file"
        and:
        def link = root.resolve('link.nf')
        Files.createSymbolicLink(link, main)
        assert !Files.isRegularFile(link, LinkOption.NOFOLLOW_LINKS)
        when:
        ResourcesBundle.scan(root)
        then:
        def e = thrown(IllegalArgumentException)
        e.message.startsWith('Module bundle does not allow link files')

    }

    def 'should find files for the given pattern' () {
        given:
        def root = folder.resolve('mod1'); root.mkdir()
        def main = root.resolve('main.nf'); main.text = "I'm the main file"
        and:
        root.resolve('bin').mkdirs()
        root.resolve('bin/hola.sh').text = 'hola'
        root.resolve('bin/sub1').mkdirs()
        root.resolve('bin/sub1/file1').text = 'file1'
        root.resolve('bin/sub2').mkdirs()
        root.resolve('bin/sub2/file2').text = 'file2'
        root.resolve('bin/sub2/file22').text = 'file22'
        and:
        root.resolve('foo').mkdirs()
        root.resolve('foo/aaa').text = 'aaa'
        and:
        root.resolve('bar').mkdirs()
        root.resolve('bar/bbb').text = 'bbb'

        when:
        def module = ResourcesBundle.scan(root, [filePattern: '{bin,bin/**}', baseDirectory: '/usr/local'])
        then:
        module.getEntries() == [
                '/usr/local/bin',
                '/usr/local/bin/hola.sh',
                '/usr/local/bin/sub1',
                '/usr/local/bin/sub1/file1',
                '/usr/local/bin/sub2',
                '/usr/local/bin/sub2/file2',
                '/usr/local/bin/sub2/file22',
        ] as Set
    }

    def 'should get bin paths' () {
        given:
        def root = folder.resolve('mod1'); root.mkdir()
        def main = root.resolve('main.nf'); main.text = "I'm the main file"
        and:
        root.resolve('bin').mkdirs()
        root.resolve('bin/hola.sh').text = 'hola'
        root.resolve('some/path').mkdirs()
        root.resolve('other/path').mkdirs()
        root.resolve('usr/bin').mkdirs()
        root.resolve('usr/local/bin').mkdirs()
        root.resolve('usr/local/foo').mkdirs()

        when:
        def module = ResourcesBundle.scan(root)
        then:
        module.content().size()>0
        and:
        module.getBinDirs() == [root.resolve('bin'),
                                root.resolve('usr/bin'),
                                root.resolve('usr/local/bin')]


    }

}
