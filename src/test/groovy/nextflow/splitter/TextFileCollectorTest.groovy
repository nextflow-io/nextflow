/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.splitter
import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Paths

import spock.lang.Specification

import static test.TestHelper.gunzip

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TextFileCollectorTest extends Specification {

    def 'test next name' () {

        given:
        def buffer = new TextFileCollector(Paths.get('.'))

        expect:
        buffer.getNextNameFor(Paths.get('/some/file.fa'),1) == Paths.get('/some/file.1.fa')
        buffer.getNextNameFor(Paths.get('/some/file.fa'),2) == Paths.get('/some/file.2.fa')
        buffer.getNextNameFor(Paths.get('/some/file.fa'),3) == Paths.get('/some/file.3.fa')
    }

    def 'test add text' () {

        given:
        def base = Files.createTempDirectory('test').resolve('sample.fasta')
        def buffer = new TextFileCollector(base)

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
        def buffer = new TextFileCollector(base)
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
        def collector = new TextFileCollector(base, Charset.defaultCharset(), true)

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


}
