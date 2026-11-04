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

package nextflow.pipeline

import java.nio.file.Files
import java.nio.file.Path

import org.eclipse.jgit.api.Git
import spock.lang.Specification
import spock.lang.TempDir

class PipelineBundleTest extends Specification {

    @TempDir
    Path tmp

    private Path write(Path dir, String name, String text='x') {
        final f = dir.resolve(name)
        Files.createDirectories(f.parent)
        f.text = text
        return f
    }

    private static List<String> names(List<Path> l) {
        l.collect { it.toString() }
    }

    def 'bundles only tracked files in a git checkout'() {
        given:
        write(tmp, 'main.nf')
        write(tmp, '.gitignore', 'work/\n')
        write(tmp, 'notes.txt')
        write(tmp, 'work/x')
        Git.init().setDirectory(tmp.toFile()).call().withCloseable { git ->
            git.add().addFilepattern('main.nf').addFilepattern('.gitignore').call()
        }

        expect:
        names(PipelineBundle.files(tmp)) == ['.gitignore', 'main.nf']
    }

    def 'honours root .gitignore outside git'() {
        given:
        write(tmp, 'main.nf')
        write(tmp, '.gitignore', 'work/\n*.log\n')
        write(tmp, 'work/a')
        write(tmp, 'x.log')

        expect:
        names(PipelineBundle.files(tmp)) == ['.gitignore', 'main.nf']
    }

    def 'skips registry-installed module dirs outside git'() {
        given:
        write(tmp, 'main.nf')
        write(tmp, 'modules/nf-core/a/main.nf')
        write(tmp, 'modules/nf-core/a/.module-info')
        write(tmp, 'modules/nf-core/b/main.nf')

        expect:
        names(PipelineBundle.files(tmp)) == ['main.nf', 'modules/nf-core/b/main.nf']
    }

    def 'never includes .pipeline-info or .git'() {
        given:
        write(tmp, 'main.nf')
        write(tmp, '.pipeline-info')

        expect:
        names(PipelineBundle.files(tmp)) == ['main.nf']

        when:
        Git.init().setDirectory(tmp.toFile()).call().withCloseable { git ->
            git.add().addFilepattern('.').call()
        }

        then:
        names(PipelineBundle.files(tmp)) == ['main.nf']
    }

    def 'create then extract round-trips the files'() {
        given:
        write(tmp.resolve('src'), 'main.nf', 'hello')
        write(tmp.resolve('src'), 'modules/local/m.nf', 'mod')
        final bundle = tmp.resolve('out/p.tar.gz')

        when:
        final sum = PipelineBundle.create(tmp.resolve('src'), bundle)
        PipelineBundle.extract(bundle, tmp.resolve('dst'))

        then:
        sum ==~ /^sha256:[0-9a-f]{64}$/
        sum == 'sha256:' + nextflow.module.ModuleStorage.computeBundleChecksum(bundle)
        tmp.resolve('dst/main.nf').text == 'hello'
        tmp.resolve('dst/modules/local/m.nf').text == 'mod'
    }

    def 'bundles a pipeline living in a repo subdirectory'() {
        given:
        write(tmp, 'other.txt')
        write(tmp, 'pipelines/p/main.nf')
        write(tmp, 'pipelines/p/untracked.txt')
        Git.init().setDirectory(tmp.toFile()).call().withCloseable { git ->
            git.add().addFilepattern('other.txt').addFilepattern('pipelines/p/main.nf').call()
        }

        expect:
        names(PipelineBundle.files(tmp.resolve('pipelines/p'))) == ['main.nf']
    }

    def 'does not bundle symlinks'() {
        given:
        write(tmp, 'main.nf')
        Files.createSymbolicLink(tmp.resolve('link.nf'), tmp.resolve('main.nf'))

        expect:
        names(PipelineBundle.files(tmp)) == ['main.nf']

        when:
        Git.init().setDirectory(tmp.toFile()).call().withCloseable { git ->
            git.add().addFilepattern('.').call()
        }

        then:
        names(PipelineBundle.files(tmp)) == ['main.nf']
    }

    def 'keeps the executable bit on round-trip'() {
        given:
        final src = tmp.resolve('src')
        write(src, 'bin/run.sh').toFile().setExecutable(true)
        write(src, 'main.nf')
        final bundle = tmp.resolve('p.tar.gz')

        when:
        PipelineBundle.create(src, bundle)
        PipelineBundle.extract(bundle, tmp.resolve('dst'))

        then:
        Files.isExecutable(tmp.resolve('dst/bin/run.sh'))
        !Files.isExecutable(tmp.resolve('dst/main.nf'))
    }

    def 'round-trips a path longer than 100 bytes'() {
        given:
        final src = tmp.resolve('src')
        final rel = 'deep/' + ('d' * 60) + '/' + ('e' * 60) + '/main.nf'
        write(src, rel, 'long')
        final bundle = tmp.resolve('p.tar.gz')

        when:
        PipelineBundle.create(src, bundle)
        PipelineBundle.extract(bundle, tmp.resolve('dst'))

        then:
        tmp.resolve('dst').resolve(rel).text == 'long'
    }

    def 'gives the same checksum for the same tree'() {
        given:
        final src = tmp.resolve('src')
        write(src, 'main.nf')
        final a = PipelineBundle.create(src, tmp.resolve('a.tar.gz'))
        sleep(1100)
        src.resolve('main.nf').toFile().setLastModified(System.currentTimeMillis())

        expect:
        a == PipelineBundle.create(src, tmp.resolve('b.tar.gz'))
    }

    def 'takes the executable bit from the git index'() {
        given:
        write(tmp, 'run.sh').toFile().setExecutable(true)
        write(tmp, 'main.nf')
        Git.init().setDirectory(tmp.toFile()).call().withCloseable { git ->
            git.add().addFilepattern('run.sh').addFilepattern('main.nf').call()
        }
        final bundle = tmp.resolve('p.tar.gz')
        tmp.resolve('run.sh').toFile().setExecutable(false)

        when:
        PipelineBundle.create(tmp, bundle)
        PipelineBundle.extract(bundle, tmp.resolve('dst'))

        then:
        Files.isExecutable(tmp.resolve('dst/run.sh'))
        !Files.isExecutable(tmp.resolve('dst/main.nf'))
    }

    def 'wraps failures in AbortOperationException and removes the partial file'() {
        given:
        final bundle = tmp.resolve('p.tar.gz')

        when:
        PipelineBundle.create(tmp.resolve('missing'), bundle)

        then:
        thrown(nextflow.exception.AbortOperationException)
        !Files.exists(bundle)
    }

    def 'extract rejects entries escaping the target'() {
        given:
        final bundle = tmp.resolve('evil.tar.gz')
        new org.apache.commons.compress.archivers.tar.TarArchiveOutputStream(
                new org.apache.commons.compress.compressors.gzip.GzipCompressorOutputStream(Files.newOutputStream(bundle))).withCloseable { tos ->
            final e = new org.apache.commons.compress.archivers.tar.TarArchiveEntry('../evil')
            e.setSize(1)
            tos.putArchiveEntry(e)
            tos.write(65)
            tos.closeArchiveEntry()
        }

        when:
        PipelineBundle.extract(bundle, tmp.resolve('dst'))

        then:
        thrown(Exception)
        !Files.exists(tmp.resolve('evil'))
    }
}
