/*
 * Copyright 2020, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ModuleBundleTest extends Specification {

    def LAST_MODIFIED = 1_000_000_000_000

    def 'should scan bundle files' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def bundlePath = folder.resolve('mod1'); bundlePath.mkdir()
        bundlePath.resolve('main.nf').text = "I'm the main file"
        bundlePath.resolve('this/that').mkdirs()
        Files.createFile(bundlePath.resolve('this/hola.txt'))
        Files.createFile(bundlePath.resolve('this/hello.txt'))
        Files.createFile(bundlePath.resolve('this/that/ciao.txt'))
        and:
        FileHelper.visitFiles([type:'any'], bundlePath, '**', { Path it -> it.setLastModified(LAST_MODIFIED) })

        when:
        def bundle = ModuleBundle.scan(bundlePath)
        then:
        bundle
        !bundle.isEmpty()
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
        bundle.fingerprint() == 'c063b8f42cfd65f4eb8efe7c30c108bc'

        cleanup:
        folder?.deleteDir()
    }

    def 'should get dockerfile' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def dockerPath = folder.resolve('Dockerfile'); dockerPath.text = "I'm the main file"
        def bundlePath = folder.resolve('bundle')
        and:
        dockerPath.setLastModified(LAST_MODIFIED)

        when:
        def bundle = ModuleBundle.scan(bundlePath)
        then:
        bundle.getDockerfile() == dockerPath
        and:
        bundle
        !bundle.isEmpty()
        and:
        bundle.fingerprint() == 'c8597047abea34f987e6a347dca36823'

        when:
        // changing file permissions, change the fingerprint
        dockerPath.setPermissions(6,0,0)
        then:
        bundle.fingerprint() == '280d52c24debce950148f4250a34e3ff'

        when:
        // changing the last modified time, change the fingerprint
        dockerPath.setLastModified(LAST_MODIFIED +100)
        then:
        bundle.fingerprint() == '41bd15592039e3a198bef861800d3cd6'
        
        cleanup:
        folder?.deleteDir()
    }

    def 'should check max file size'() {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def root = folder.resolve('mod1'); root.mkdir()
        root.resolve('main.nf').text = "I'm the main file"

        when:
        ModuleBundle.scan(root, [maxFileSize: MemoryUnit.of(5)])
        then:
        thrown(IllegalArgumentException)

        cleanup:
        folder?.deleteDir()
    }

    def 'should check max bundle size'() {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def root = folder.resolve('mod1'); root.mkdir()
        root.resolve('main.nf').text = "I'm the main file"

        when:
        ModuleBundle.scan(root, [maxBundleSize: MemoryUnit.of(5)])
        then:
        def e = thrown(IllegalArgumentException)
        e.message == 'Module total size cannot exceed 5 B'

        cleanup:
        folder?.deleteDir()
    }

    def 'should symlink not allowd'() {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def root = folder.resolve('mod1'); root.mkdir()
        def main = root.resolve('main.nf'); main.text = "I'm the main file"
        and:
        def link = root.resolve('link.nf')
        Files.createSymbolicLink(link, main)
        assert !Files.isRegularFile(link, LinkOption.NOFOLLOW_LINKS)
        when:
        ModuleBundle.scan(root)
        then:
        def e = thrown(IllegalArgumentException)
        e.message.startsWith('Module bundle does not allow link files')

        cleanup:
        folder?.deleteDir()
    }
}
