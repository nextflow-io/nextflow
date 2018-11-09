/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.splitter
import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Paths

import spock.lang.Specification

import com.google.common.hash.HashCode
import static test.TestHelper.gunzip

import nextflow.splitter.TextFileCollector.CachePath

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TextFileCollectorTest extends Specification {

    def 'test next name' () {

        given:
        def base = new TextFileCollector.CachePath(Paths.get('.'))
        def buffer = new TextFileCollector(base)

        expect:
        buffer.getNextNameFor(Paths.get('/some/file.fa'),1) == Paths.get('/some/file.1.fa')
        buffer.getNextNameFor(Paths.get('/some/file.fa'),2) == Paths.get('/some/file.2.fa')
        buffer.getNextNameFor(Paths.get('/some/file.fa'),3) == Paths.get('/some/file.3.fa')
    }

    def 'test add text' () {

        given:
        def base = Files.createTempDirectory('test').resolve('sample.fasta')
        def buffer = new TextFileCollector(new CachePath(base))

        when:
        buffer.add('>seq1\n')
        buffer.add('alpha\n')
        assert buffer.nextChunk() == base.resolveSibling('sample.1.fasta')

        buffer.add('>seq2\n')
        buffer.add('gamma\n')
        buffer.add('>seq3\n')
        buffer.add('beta\n')
        assert buffer.nextChunk() == base.resolveSibling('sample.2.fasta')

        buffer.add('>seq4\n')
        buffer.add('kappa\n')
        buffer.add('>seq5\n')
        buffer.add('iota\n')
        buffer.add('delta\n')
        assert buffer.nextChunk() == base.resolveSibling('sample.3.fasta')

        buffer.close()

        then:
        base.resolveSibling('sample.1.fasta').text == '>seq1\nalpha\n'
        base.resolveSibling('sample.2.fasta').text == '>seq2\ngamma\n>seq3\nbeta\n'
        base.resolveSibling('sample.3.fasta').text == '>seq4\nkappa\n>seq5\niota\ndelta\n'
        base.resolveSibling('.chunks.sample.fasta').exists()
        buffer.checkCached()
        buffer.getAllChunks()*.name == ['sample.1.fasta','sample.2.fasta','sample.3.fasta']

        cleanup:
        base?.parent?.deleteDir()
    }


    def 'test isEmpty' () {

        given:
        def base = Files.createTempDirectory('test').resolve('chunk.fasta')

        when:
        def buffer = new TextFileCollector(new CachePath(base))
        then:
        !buffer.hasChunk()

        when:
        buffer.add('>seq1\n')
        buffer.add('alpha\n')
        then:
        buffer.hasChunk()

        when:
        buffer.nextChunk()
        buffer.add('>seq2\n')
        buffer.add('gamma\n')
        buffer.add('>seq3\n')
        buffer.add('beta\n')
        then:
        buffer.hasChunk()

        when:
        buffer.nextChunk()
        then:
        !buffer.hasChunk()

        cleanup:
        base?.parent?.deleteDir()
    }

    def 'test compress chunks' () {

        given:
        def base = Files.createTempDirectory('test').resolve('sample.fasta')
        def collector = new TextFileCollector(new CachePath(base), Charset.defaultCharset(), true)

        when:
        collector.add('>seq1\n')
        collector.add('alpha\n')
        assert collector.nextChunk() == base.resolveSibling('sample.1.fasta.gz')

        collector.add('>seq2\n')
        collector.add('gamma\n')
        collector.add('>seq3\n')
        collector.add('beta\n')
        assert collector.nextChunk() == base.resolveSibling('sample.2.fasta.gz')

        collector.add('>seq4\n')
        collector.add('kappa\n')
        collector.add('>seq5\n')
        collector.add('iota\n')
        collector.add('delta\n')
        assert collector.nextChunk() == base.resolveSibling('sample.3.fasta.gz')

        collector.close()

        then:
        gunzip(base.resolveSibling('sample.1.fasta.gz')) == '>seq1\nalpha\n'
        gunzip(base.resolveSibling('sample.2.fasta.gz')) == '>seq2\ngamma\n>seq3\nbeta\n'
        gunzip(base.resolveSibling('sample.3.fasta.gz')) == '>seq4\nkappa\n>seq5\niota\ndelta\n'
        base.resolveSibling('.chunks.sample.fasta').exists()
        collector.checkCached()
        collector.getAllChunks()*.name == ['sample.1.fasta.gz','sample.2.fasta.gz','sample.3.fasta.gz']

        cleanup:
        base?.parent?.deleteDir()
    }

    def 'should create unique cache marker file' () {
        given:
        def hash = HashCode.fromLong(89323923l)
        def base = Files.createTempDirectory('test').resolve('sample.fasta')
        def collector = new TextFileCollector(new CachePath(base, hash), Charset.defaultCharset(), false)

        when:
        collector.add('xxx')
        assert collector.nextChunk() == base.resolveSibling('sample.1.fasta')

        collector.add('yyy')
        assert collector.nextChunk() == base.resolveSibling('sample.2.fasta')

        collector.add('zzz')
        assert collector.nextChunk() == base.resolveSibling('sample.3.fasta')

        collector.close()

        then:
        base.resolveSibling(".chunks.sample.fasta.$hash").exists()
        collector.checkCached()
        collector.getAllChunks()*.name == ['sample.1.fasta','sample.2.fasta','sample.3.fasta']

        cleanup:
        base?.parent?.deleteDir()
    }

}
