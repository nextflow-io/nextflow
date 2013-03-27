/*
 * Copyright (c) 2012, the authors.
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

package util

import nextflow.Nextflow
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ChunkMethodsTest extends Specification {

    def 'test chunk string by line' () {

        when:
        def out = "Hello".chunkLines()

        then:
        Nextflow.read(out)  == 'Hello'


        when:
        out = "Hello\nHola\nHalo".chunkLines()

        then:
        Nextflow.read(out) == 'Hello'
        Nextflow.read(out) == 'Hola'
        Nextflow.read(out) == 'Halo'


        when:
        out = "11\n22\n33\n44\n55".chunkLines(3)

        then:
        Nextflow.read(out) == '11\n22\n33'
        Nextflow.read(out) == '44\n55'

    }

    def 'test chunk file by line ' () {

        setup:
        def file = File.createTempFile('chunk','test')
        file.deleteOnExit()
        file.text = '''\
        line1
        line2
        line3
        line4
        line5
        '''.stripIndent()

        when:
        def channel = file.chunkLines()

        then:
        Nextflow.read(channel) == 'line1'
        Nextflow.read(channel) == 'line2'
        Nextflow.read(channel) == 'line3'
        Nextflow.read(channel) == 'line4'
        Nextflow.read(channel) == 'line5'

        when:
        channel = file.chunkLines(2)

        then:
        Nextflow.read(channel) == 'line1\nline2'
        Nextflow.read(channel) == 'line3\nline4'
        Nextflow.read(channel) == 'line5'

    }


    def 'test check by Fasta ' () {

        when:
        def fasta = """\
                >prot1
                LCLYTHIGRNIYYGS1
                EWIWGGFSVDKATLN
                >prot2
                LLILILLLLLLALLS
                GLMPFLHTSKHRSMM
                IENY
                """.stripIndent()

        def result = fasta.chunkFasta()

        then:
        Nextflow.read(result) == ">prot1\nLCLYTHIGRNIYYGS1\nEWIWGGFSVDKATLN\n"
        Nextflow.read(result) == ">prot2\nLLILILLLLLLALLS\nGLMPFLHTSKHRSMM\nIENY\n"

    }

    def 'test check by Fasta file' () {

        setup:
        def file = File.createTempFile('chunk','test')
        file.deleteOnExit()
        def fasta = """\
                >prot1
                AA
                >prot2
                BB
                CC
                >prot3
                DD
                >prot4
                EE
                FF
                GG
                >prot5
                LL
                NN
                """.stripIndent()


        when:
        def result = fasta.chunkFasta(2)

        then:
        Nextflow.read(result) == ">prot1\nAA\n>prot2\nBB\nCC\n"
        Nextflow.read(result) == ">prot3\nDD\n>prot4\nEE\nFF\nGG\n"
        Nextflow.read(result) == ">prot5\nLL\nNN\n"

    }

}
