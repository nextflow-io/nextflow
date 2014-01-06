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

package nextflow.extension

import java.nio.charset.Charset

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.operator.PoisonPill
import groovyx.gpars.group.NonDaemonPGroup
import nextflow.Nextflow
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


    def testChopInvoke() {

        expect:
        NextflowExtensions.chopInvoke(null, 'hola', 1) == 'hola'
        NextflowExtensions.chopInvoke({ x -> x.reverse() }, 'hola', 2) == 'aloh'
        NextflowExtensions.chopInvoke({ x, y -> y * 2 }, 'hola', 3) == 6

    }

    def testFastaRecord() {
        def fasta = /
            ;
            >1aboA  xyz|beta|gamma

            NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS
            NYITPVN
            /.stripIndent()

        expect:
        NextflowExtensions.parseFastaRecord(fasta, [id:true]) .id == '1aboA'
        NextflowExtensions.parseFastaRecord(fasta, [seq:true]) .seq == 'NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS\nNYITPVN'
        NextflowExtensions.parseFastaRecord(fasta, [seq:true, width: 20 ]) .seq == 'NLFVALYDFVASGDNTLSIT\nKGEKLRVLGYNHNGEWCEAQ\nTKNGQGWVPSNYITPVN'
        NextflowExtensions.parseFastaRecord(fasta, [head:true]) .head == '1aboA  xyz|beta|gamma'
        NextflowExtensions.parseFastaRecord(fasta, [desc:true, seq: true, width: 0]) .seq == 'NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVN'

    }


    def 'test chop text by line' () {

        expect:
        // no params are specified, it returns the last line chopped
        "Hello".chopLines() == 'Hello'
        "Hello\nworld".chopLines() == 'world'
        // when *into* is a List object, it return that list populated
        "Hello\nworld\n!".chopLines(into:[]) == ['Hello','world','!']
        "Hello\nworld\n!".chopLines(into:[]) { it.reverse() } == ['olleH','dlrow','!']

        when:   '*into* is a Dataflow channel, get the chopped items, closing by a poison-pill'
        def channel = "Hello\nworld\n!".chopLines(into: new DataflowQueue())
        then:
        channel.val == 'Hello'
        channel.val == 'world'
        channel.val == '!'
        channel.val == PoisonPill.instance


        expect:
        "Hello\nHola\nHalo".chopLines(into:[]) == ['Hello', 'Hola', 'Halo']
        "11\n22\n33\n44\n55".chopLines(count:3, into:[]) == [ '11\n22\n33', '44\n55' ]
        "11\n22\n33\n44\n55".chopLines(count:2, into:[]) == [ '11\n22', '33\n44', '55' ]

    }

    def 'test chop file by line' () {

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
        def lines = file.chopLines(into:[])

        then:
        lines[0] == 'line1'
        lines[1] == 'line2'
        lines[2]== 'line3'
        lines[3] == 'line4'
        lines[4] == 'line5'

        when:
        def channel = new DataflowQueue()
        file.chopLines(count:2, into: channel)

        then:
        channel.val == 'line1\nline2'
        channel.val == 'line3\nline4'
        channel.val == 'line5'

    }

    def 'test chop Fasta ' () {

        when:
        def fasta = """\
                >prot1
                LCLYTHIGRNIYYGS1
                EWIWGGFSVDKATLN
                ;
                ; comment
                ;
                >prot2
                LLILILLLLLLALLS
                GLMPFLHTSKHRSMM
                IENY
                """.stripIndent()

        def count = 0
        def q = fasta.chopFasta(into: new DataflowQueue()) { count++; it }

        then:
        count == 2
        q.val == ">prot1\nLCLYTHIGRNIYYGS1\nEWIWGGFSVDKATLN\n"
        q.val == ">prot2\nLLILILLLLLLALLS\nGLMPFLHTSKHRSMM\nIENY\n"
        q.val == PoisonPill.instance

    }

    def 'test chop Fasta record ' () {

        when:
        def fasta = """\
                >1aboA
                NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS
                NYITPVN
                >1ycsB
                KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY
                VPRNLLGLYP
                ; comment
                >1pht
                GYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG
                WLNGYNETTGERGDFPGTYVE
                YIGRKKISP
                """.stripIndent()

        def q = fasta.chopFasta(into: new DataflowQueue(), record: [id:true, seq:true, width: 999])

        then:
        q.val == [id:'1aboA', seq: 'NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVN']
        q.val == [id:'1ycsB', seq: 'KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGYVPRNLLGLYP']
        q.val == [id:'1pht', seq: 'GYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIGWLNGYNETTGERGDFPGTYVEYIGRKKISP']
        q.val == PoisonPill.instance

    }

    def 'test chop Fasta file' () {

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
        def result = fasta.chopFasta(count:2, into: [])

        then:
        result[0] == ">prot1\nAA\n>prot2\nBB\nCC\n"
        result[1] == ">prot3\nDD\n>prot4\nEE\nFF\nGG\n"
        result[2] == ">prot5\nLL\nNN\n"

    }

    def 'test chop string' () {

        expect:
        '012345678901234567'.chopString(count:5, into:[]) == ['01234','56789','01234','567']
        '012345678901234567'.chopString(count:5, into:[]) {it.reverse()} == ['43210','98765','43210','765']

        when:
        def q = '012345678901234567'.chopString(count:5, into:new DataflowQueue())
        then:
        q.val == '01234'
        q.val == '56789'
        q.val == '01234'
        q.val == '567'
        q.val == PoisonPill.instance
    }


    def 'test chop bytes' () {

        setup:
        def bytes = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6 ] as byte[]

        expect:
        new ByteArrayInputStream(bytes).chopBytes(count: 10) == [0, 1, 2, 3, 4, 5, 6] as byte[]
        new ByteArrayInputStream(bytes).chopBytes(count:5, into:[]) == [ [0, 1, 2, 3, 4] as byte[], [5, 6, 7, 8, 9] as byte[] , [ 0, 1, 2, 3, 4] as byte[], [5, 6] as byte[] ]

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

    def testGetCharSet() {

        expect:
        NextflowExtensions.getCharset('x') == Charset.defaultCharset()
        NextflowExtensions.getCharset([x:1]) == Charset.defaultCharset()
        NextflowExtensions.getCharset('iso-8859-1') == Charset.forName('iso-8859-1')
        NextflowExtensions.getCharset(charset:'iso-8859-1') == Charset.forName('iso-8859-1')
        NextflowExtensions.getCharset(xx:'iso-8859-1') == Charset.defaultCharset()

    }

}
