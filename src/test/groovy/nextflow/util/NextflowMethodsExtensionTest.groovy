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

package nextflow.util

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.Nextflow
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NextflowMethodsExtensionTest extends Specification {

    def 'test chunk string by line' () {

        when:
        def out = []
        "Hello".chunkLines { out << it }

        then:
        out  == ['Hello']


        when:
        out = []
        "Hello\nHola\nHalo".chunkLines { out << it }

        then:
        out == ['Hello', 'Hola', 'Halo']


        when:
        out = []
        "11\n22\n33\n44\n55".chunkLines(3) { out << it }

        then:
        out == [ '11\n22\n33', '44\n55' ]

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
        def channel = new DataflowQueue()
        file.chunkLines { channel << it }

        then:
        Nextflow.read(channel) == 'line1'
        Nextflow.read(channel) == 'line2'
        Nextflow.read(channel) == 'line3'
        Nextflow.read(channel) == 'line4'
        Nextflow.read(channel) == 'line5'

        when:
        channel = new DataflowQueue()
        file.chunkLines(2) { channel << it }

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

        def result = new DataflowQueue()
        fasta.chunkFasta { result << it }

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
        def result = new DataflowQueue()
        fasta.chunkFasta(2) { result << it }

        then:
        Nextflow.read(result) == ">prot1\nAA\n>prot2\nBB\nCC\n"
        Nextflow.read(result) == ">prot3\nDD\n>prot4\nEE\nFF\nGG\n"
        Nextflow.read(result) == ">prot5\nLL\nNN\n"

    }


    def data1 = new DataflowQueue<>()
    def data2 = new DataflowQueue()

    def 'test closure interceptor' () {

        setup:

        def closure = {

            data1 << 1
            data2

        }

        when:
        def interceptor = new NextflowMethodsExtension.WritableChannelInterceptor(closure.owner)
        closure.@owner = interceptor
        closure.call()

        then:
        interceptor.getWrittenChannels() *.unwrap() == [data1]

    }


    def channel1 = new DataflowQueue<>()
    def channel2 = new DataflowQueue()

    def 'test Each' () {

        setup:
        def queue = new DataflowQueue()

        def closure = { it ->
            channel1 << it
            channel2 << it * it

        }

        when:
        NextflowMethodsExtension.each(queue, closure);
        queue << 1 << 2 << PoisonPill

        then:
        channel1.getVal() == 1
        channel1.getVal() == 2
        channel2.getVal() == 1
        channel2.getVal() == 4

    }

    def 'test cloure ' () {

        when:
        def closure = {}
        println "delegate: ${closure.delegate}"
        println "owner: ${closure.owner}"


        then: true
    }

}
