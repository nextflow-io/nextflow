/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.extension
import java.nio.file.Path
import java.nio.file.Paths

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.group.NonDaemonPGroup
import nextflow.util.Duration
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NextflowExtensionsTest extends Specification {


    def 'test leftTrim' () {

        expect:
        '  hola hello  '.leftTrim() == 'hola hello  '
        '\n\n hola hello\n'.leftTrim() == 'hola hello\n'

    }

    def 'test rightTrim' () {

        expect:
        '  hola hello  '.rightTrim() == '  hola hello'
        '\n\nhola hello\n\n'.rightTrim() == '\n\nhola hello'

    }



    def 'test chunk string by line' () {

        when:
        def out = []; "Hello".chunkLines { out << it }
        then:
        out  == ['Hello']

        when:
        out = []; "Hello\nHola\nHalo".chunkLines { out << it }
        then:
        out == ['Hello', 'Hola', 'Halo']

        when:
        out = []; "11\n22\n33\n44\n55".chunkLines(3) { out << it }
        then:
        out == [ '11\n22\n33', '44\n55' ]

        when:
        out = []; "11\n22\n33\n44\n55".chunkLines(size: 2) { out << it }
        then:
        out == [ '11\n22', '33\n44', '55' ]

    }

    // the channel has to be defined at class level, otherwise the 'WritableChannelInterceptor' does not work
    def channel = new DataflowQueue()

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
        file.chunkLines { channel << it }

        then:
        channel.val == 'line1'
        channel.val == 'line2'
        channel.val == 'line3'
        channel.val == 'line4'
        channel.val == 'line5'


        when:
        channel = new DataflowQueue()
        def closure = { channel << it }
        file.chunkLines(2, closure)

        then:
        channel.val == 'line1\nline2'
        channel.val == 'line3\nline4'
        channel.val == 'line5'
        // the following getter are defined by the 'WritableChannelInterceptor'
        closure.@owner.getWrittenChannels().size() == 1
        closure.@owner.isClosed() == true


        when:
        channel = new DataflowQueue()
        closure = { channel << it }
        file.chunkLines(size:3, autoClose: false, closure)

        then:
        channel.val == 'line1\nline2\nline3'
        channel.val == 'line4\nline5'
        // the following getter are defined by the 'WritableChannelInterceptor'
        closure.@owner.getWrittenChannels().size() == 1 // one channel has been intercepted
        closure.@owner.isClosed() == false              // no channel has been close as requested


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
        result.val == ">prot1\nLCLYTHIGRNIYYGS1\nEWIWGGFSVDKATLN\n"
        result.val == ">prot2\nLLILILLLLLLALLS\nGLMPFLHTSKHRSMM\nIENY\n"

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
        result.val == ">prot1\nAA\n>prot2\nBB\nCC\n"
        result.val == ">prot3\nDD\n>prot4\nEE\nFF\nGG\n"
        result.val == ">prot5\nLL\nNN\n"

    }


    def testPdbSeqRef() {

        setup:
        def pdb = /
        >1bvh_A mol:protein length:157  ACID PHOSPHATASE
        AEQVTKSVLFVCLGNICRSPIAEAVFRKLVTDQNISDNWVIDSGAVSDWNVGRSPDPRAVSCLRNHGINTAHKARQVTKEDFVTFDYILCMDESNLRDLNRKSNQVKNCRAKIELLGSYDPQKQLIIEDPYYGNDADFETVYQQCVRCCR
        AFLEKVR
        >1bvi_A mol:protein length:104  PROTEIN (RIBONUCLEASE T1)
        ACDYTCGSNCYSSSDVSTAQAAGYKLHEDGETVGSNSYPHKYNNYEGFDFSVSSPYYEWPILSSGDVYSGGSPGADRVVFNENNQLAGVITHTGASGNNFVECT
        >1bvi_B mol:protein length:104  PROTEIN (RIBONUCLEASE T1)
        ACDYTCGSNCYSSSDVSTAQAAGYKLHEDGETVGSNSYPHKYNNYEGFDFSVSSPYYEWPILSSGDVYSGGSPGADRVVFNENNQLAGVITHTGASGNNFVECT
        >1bvi_C mol:protein length:104  PROTEIN (RIBONUCLEASE T1)
        ACDYTCGSNCYSSSDVSTAQAAGYKLHEDGETVGSNSYPHKYNNYEGFDFSVSSPYYEWPILSSGDVYSGGSPGADRVVFNENNQLAGVITHTGASGNNFVECT
        >1bvi_D mol:protein length:104  PROTEIN (RIBONUCLEASE T1)
        ACDYTCGSNCYSSSDVSTAQAAGYKLHEDGETVGSNSYPHKYNNYEGFDFSVSSPYYEWPILSSGDVYSGGSPGADRVVFNENNQLAGVITHTGASGNNFVECT
        >1bvk_A mol:protein length:108  HULYS11
        DIQMTQSPSSLSASVGDRVTITCRASGNIHNYLAWYQQKPGKAPKLLIYYTTTLADGVPSRFSGSGSGTDYTFTISSLQPEDIATYYCQHFWSTPRTFGQGTKVEIKR
        /.stripIndent()

        when:
        def result = (List)NextflowExtensions.chopFasta(pdb, [into:[], record: [id:true] ])
        then:
        result.size() == 6
        result[0].id == '1bvh_A'

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
        def interceptor = new NextflowExtensions.WritableChannelInterceptor(closure)
        closure.call()

        then:
        interceptor.getWrittenChannels() *.unwrap() == [data1]

    }


    def channel1 = new DataflowQueue()
    def channel2 = new DataflowQueue()

    def 'test Each' () {

        setup:
        def queue = new DataflowQueue()

        def closure = { it ->
            channel1 << it
            channel2 << it * it
        }

        when:
        // note: launching the test from Gradle, it requires the NonDaemonPGroup to be specified, otherwise it raises a RejectedExecutionException exception
        NextflowExtensions.each(queue, new NonDaemonPGroup(), closure);
        queue << 1 << 2

        then:
        channel1.getVal() == 1
        channel1.getVal() == 2
        channel2.getVal() == 1
        channel2.getVal() == 4

    }

    def 'test closure ' () {

        when:
        def closure = {}
        println "delegate: ${closure.delegate}"
        println "owner: ${closure.owner}"


        then: true
    }

    def testAsDuration() {

        setup:
        def x = 3;

        expect:
        2_000 as Duration == Duration.of('2 second')
        '1s' as Duration == Duration.of('1 second')
        "$x min" as Duration == Duration.of('3 min')

    }

    def testAsPath() {

        setup:
        def x = 'Hello'

        expect:
        'file.txt' as Path == Paths.get('file.txt')
        '/some/path/file.txt' as Path == Paths.get('/some/path/file.txt')
        "name.fa" as Path == Paths.get('name.fa')
        "/some/path/${x}.txt" as Path == Paths.get('/some/path/Hello.txt')

        new File('/path/to/file') as Path == Paths.get('/path/to/file')


    }



}
