/*
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
 */

package nextflow.splitter

import spock.lang.Specification

import java.nio.file.Files
import java.util.zip.GZIPOutputStream

import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SplitterExTest extends Specification {


    def 'should count lines'() {

        given:
        def file = TestHelper.createInMemTempFile('foo.txt')

        when:
        file.text = '''
            line 1
            line 2
            line 3
            line 4
            line 5
            '''
                .stripIndent().strip()

        then:
        file.countText() == 5
        file.countText(by:2) == 3
    }

    def 'should count compressed text file' () {
        given:
        def folder = TestHelper.createInMemTempDir()
        def file = folder.resolve('file.txt.gz')
        def out = new GZIPOutputStream(Files.newOutputStream(file))
        out << '''
            line 1
            line 2
            line 3
            line 4
            line 5
            line 6 
            line 7
            line 8
            '''
                .stripIndent().strip()
        out.close()

        expect:
        file.countLines() == 8
    }

    def 'should count fasta records'() {

        given:
        def file = TestHelper.createInMemTempFile('foo.txt')
        file.text = '''
            >1aboA
            NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS
            NYITPVN
            >1ycsB
            KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY
            VPRNLLGLYP
            >1pht
            GYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG
            WLNGYNETTGERGDFPGTYVEYIGRKKISP
            >1vie
            DRVRKKSGAAWQGQIVGWYCTNLTPEGYAVESEAHPGSVQIYPVAALERI
            N
            >1ihvA
            NFRVYYRDSRDPVWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRD
            '''
                .stripIndent()

        expect:
        file.countFasta() == 5
        file.countFasta(by:4) == 2
    }

    def 'should count compressed fasta' () {
        given:
        def folder = TestHelper.createInMemTempDir()
        def file = folder.resolve('fasta.gz')
        def out = new GZIPOutputStream(Files.newOutputStream(file))
        out << '''
            >1aboA
            NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS
            NYITPVN
            >1ycsB
            KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY
            VPRNLLGLYP
            >1pht
            GYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG
            WLNGYNETTGERGDFPGTYVEYIGRKKISP
            >1vie
            DRVRKKSGAAWQGQIVGWYCTNLTPEGYAVESEAHPGSVQIYPVAALERI
            N
            '''.stripIndent().strip()
        out.close()

        expect:
        file.countFasta() == 4

    }

    def 'should count fastq records' () {
        given:
        def file = TestHelper.createInMemTempFile('foo.txt')
        file.text = '''
        @SRR636272.19519409/1
        GGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGATCG
        +
        CCCFFFFDHHD;FF=GGDHGGHIIIGHIIIBDGBFCAHG@E=6?CBDBB;?BB@BD8BB;BDB<>>;@?BB<9>&5<?288AAABDBBBBACBCAC?@AD?CAC?
        @SRR636272.13995011/1
        GCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAGATCGGAAGAGCACACGTCTGAACTCC
        +
        BBCFDFDEFFHHFIJIHGHGHGIIFIJJJJIGGBFHHIEGBEFEFFCDDDD:@@<BB8BBDDDDDDBBB?AA?CDABDD5?CDDDBB<A<>ACBB8ACDCD@CD>
        @SRR636272.21107783/1
        CGGGGAGCGCGGGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTAGGGAGCACCCCGCCGCAGGGGGACAGGCGAGATCGGAAGAGCACACGTCT
        +
        BCCFFDFFHHHHHJJJJJIJHHHHFFFFEEEEEEEDDDDDDBDBDBBDBBDBBB(:ABCDDDDDDDDDDDDDDDD@BBBDDDDDDDDDDDDBDDDDDDDDDDADC
        @SRR636272.23331539/1
        GGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGCTGGAGGAGCTGCTGGCCGGGAGGGACTTCACCGGCGAGATCGGAAGAG
        +
        CCCFFFFFHHHHHJJJJJJJJJJJJJJJHFDDBDDBDDDDDDDDDDDDADDDDDDDDDDDDDDDDDDDDDDDDDDBDBDDD9@DDDDDDDDDDDDBBDDDBDD@@
        '''.stripIndent().leftTrim()

        expect:
        file.countFastq() == 4
        file.countFastq(by:2) == 2
    }

    def 'should count compressed fastq file' () {
        given:
        def folder = TestHelper.createInMemTempDir()
        def file = folder.resolve('fastq.gz')
        def out = new GZIPOutputStream(Files.newOutputStream(file))
        out << '''
            @SRR636272.19519409/1
            GGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGATCG
            +
            CCCFFFFDHHD;FF=GGDHGGHIIIGHIIIBDGBFCAHG@E=6?CBDBB;?BB@BD8BB;BDB<>>;@?BB<9>&5<?288AAABDBBBBACBCAC?@AD?CAC?
            @SRR636272.13995011/1
            GCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAGATCGGAAGAGCACACGTCTGAACTCC
            +
            BBCFDFDEFFHHFIJIHGHGHGIIFIJJJJIGGBFHHIEGBEFEFFCDDDD:@@<BB8BBDDDDDDBBB?AA?CDABDD5?CDDDBB<A<>ACBB8ACDCD@CD>
            @SRR636272.21107783/1
            CGGGGAGCGCGGGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTAGGGAGCACCCCGCCGCAGGGGGACAGGCGAGATCGGAAGAGCACACGTCT
            +
            BCCFFDFFHHHHHJJJJJIJHHHHFFFFEEEEEEEDDDDDDBDBDBBDBBDBBB(:ABCDDDDDDDDDDDDDDDD@BBBDDDDDDDDDDDDBDDDDDDDDDDADC
            @SRR636272.23331539/1
            GGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGCTGGAGGAGCTGCTGGCCGGGAGGGACTTCACCGGCGAGATCGGAAGAG
            +
            CCCFFFFFHHHHHJJJJJJJJJJJJJJJHFDDBDDBDDDDDDDDDDDDADDDDDDDDDDDDDDDDDDDDDDDDDDBDBDDD9@DDDDDDDDDDDDBBDDDBDD@@
            '''.stripIndent().strip()
        out.close()
        
        expect:
        file.countFastq() == 4

    }

    def 'should split text file'() {
        given:
        def file = TestHelper.createInMemTempFile('foo.txt')
        file.text = "hello\nworld"

        when:
        def list = file.splitText()
        then:
        list == ['hello\n','world\n']

    }

    def 'should split text and invoke each'() {
        given:
        def file = TestHelper.createInMemTempFile('foo.txt')
        file.text = "hello\nworld"
        
        when:
        def list = file.splitText(each: { it.trim().reverse() })
        then:
        list == ['olleh','dlrow']

    }


    def 'should split text with options'() {

        given:
        def file = TestHelper.createInMemTempFile('foo.txt')

        when:
        file.text = 'Hello\nworld\n!'
        def list = file.splitText()
        then:
        list == ['Hello\n','world\n','!\n']

        when:
        file.text = 'Hello\nworld\n!'
        list = file.splitText( by:2)
        then:
        list == ['Hello\nworld\n','!\n']

        when:
        file.text = 'Hello\nworld\n!'
        list = file.splitText(each: { it.trim().reverse() })
        then:
        list == ['olleH','dlrow','!']

    }


    def 'should split fasta records' () {

        given:
        def file = TestHelper.createInMemTempFile('foo.txt')
        file.text = '''\
            >1aboA
            NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS
            NYITPVN
            >1ycsB
            KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY
            VPRNLLGLYP
            >1pht
            GYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG
            WLNGYNETTGERGDFPGTYVEYIGRKKISP
            >1vie
            DRVRKKSGAAWQGQIVGWYCTNLTPEGYAVESEAHPGSVQIYPVAALERI
            N
            >1ihvA
            NFRVYYRDSRDPVWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRD
            '''
                .stripIndent()

        when:
        def result = file.splitFasta()
        then:
        result.size() == 5 
        result[0] == '''\
            >1aboA
            NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS
            NYITPVN
            '''.stripIndent()

        result[1] == '''\
            >1ycsB
            KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY
            VPRNLLGLYP
            '''.stripIndent()

        result[2] == '''\
            >1pht
            GYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG
            WLNGYNETTGERGDFPGTYVEYIGRKKISP
            '''.stripIndent()

        result[3] == '''\
            >1vie
            DRVRKKSGAAWQGQIVGWYCTNLTPEGYAVESEAHPGSVQIYPVAALERI
            N
            '''.stripIndent()

        result[4] == '''\
            >1ihvA
            NFRVYYRDSRDPVWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRD
            '''.stripIndent()

    }


    def 'should split fastq records' () {
        given:
        def file = TestHelper.createInMemTempFile('foo.txt')
        file.text = '''\
        @SRR636272.19519409/1
        GGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGATCG
        +
        CCCFFFFDHHD;FF=GGDHGGHIIIGHIIIBDGBFCAHG@E=6?CBDBB;?BB@BD8BB;BDB<>>;@?BB<9>&5<?288AAABDBBBBACBCAC?@AD?CAC?
        @SRR636272.13995011/1
        GCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAGATCGGAAGAGCACACGTCTGAACTCC
        +
        BBCFDFDEFFHHFIJIHGHGHGIIFIJJJJIGGBFHHIEGBEFEFFCDDDD:@@<BB8BBDDDDDDBBB?AA?CDABDD5?CDDDBB<A<>ACBB8ACDCD@CD>
        @SRR636272.21107783/1
        CGGGGAGCGCGGGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTAGGGAGCACCCCGCCGCAGGGGGACAGGCGAGATCGGAAGAGCACACGTCT
        +
        BCCFFDFFHHHHHJJJJJIJHHHHFFFFEEEEEEEDDDDDDBDBDBBDBBDBBB(:ABCDDDDDDDDDDDDDDDD@BBBDDDDDDDDDDDDBDDDDDDDDDDADC
        @SRR636272.23331539/1
        GGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGCTGGAGGAGCTGCTGGCCGGGAGGGACTTCACCGGCGAGATCGGAAGAG
        +
        CCCFFFFFHHHHHJJJJJJJJJJJJJJJHFDDBDDBDDDDDDDDDDDDADDDDDDDDDDDDDDDDDDDDDDDDDDBDBDDD9@DDDDDDDDDDDDBBDDDBDD@@
        '''.stripIndent().leftTrim()

        when:
        def records = file.splitFastq()
        then:
        records.size() == 4
        records[0] == '''\
        @SRR636272.19519409/1
        GGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGATCG
        +
        CCCFFFFDHHD;FF=GGDHGGHIIIGHIIIBDGBFCAHG@E=6?CBDBB;?BB@BD8BB;BDB<>>;@?BB<9>&5<?288AAABDBBBBACBCAC?@AD?CAC?
        '''.stripIndent()

        records[1] == '''\
        @SRR636272.13995011/1
        GCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTGGGGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAGATCGGAAGAGCACACGTCTGAACTCC
        +
        BBCFDFDEFFHHFIJIHGHGHGIIFIJJJJIGGBFHHIEGBEFEFFCDDDD:@@<BB8BBDDDDDDBBB?AA?CDABDD5?CDDDBB<A<>ACBB8ACDCD@CD>
        '''.stripIndent()

        records[2] == '''\
        @SRR636272.21107783/1
        CGGGGAGCGCGGGCCCGGCAGCAGGATGATGCTCTCCCGGGCCAAGCCGGCTGTAGGGAGCACCCCGCCGCAGGGGGACAGGCGAGATCGGAAGAGCACACGTCT
        +
        BCCFFDFFHHHHHJJJJJIJHHHHFFFFEEEEEEEDDDDDDBDBDBBDBBDBBB(:ABCDDDDDDDDDDDDDDDD@BBBDDDDDDDDDDDDBDDDDDDDDDDADC
        '''.stripIndent()

        records[3] == '''\
        @SRR636272.23331539/1
        GGAGCACCCCGCCGCAGGGGGACAGGCGGAGGAAGAAAGGGAAGAAGGTGCCACAGCTGGAGGAGCTGCTGGCCGGGAGGGACTTCACCGGCGAGATCGGAAGAG
        +
        CCCFFFFFHHHHHJJJJJJJJJJJJJJJHFDDBDDBDDDDDDDDDDDDADDDDDDDDDDDDDDDDDDDDDDDDDDBDBDDD9@DDDDDDDDDDDDBBDDDBDD@@
        '''.stripIndent()
    }


    def 'should split csv' () {
        given:
        def file = TestHelper.createInMemTempFile('foo.txt')
        file.text = '''\
        a,b,c
        p,q,s
        x,y,z
        '''.stripIndent()

        when:
        def result = file.splitCsv()
        then:
        result.size() == 3 
        result[0] == ['a','b','c']
        result[1] == ['p','q','s']
        result[2] == ['x','y','z']
    }

}
