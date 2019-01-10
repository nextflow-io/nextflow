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

import nextflow.Channel
import spock.lang.Specification
import test.TestHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FastaSplitterTest extends Specification {


    def testFastaRecord() {
        def fasta = /
            ;

            >gi|5524211|gb|AAD44166.1| cytochrome b
            NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS
            NYITPVN
            /.stripIndent()

        expect:
        FastaSplitter.parseFastaRecord(fasta, [id:true])
                .id == 'gi|5524211|gb|AAD44166.1|'

        FastaSplitter.parseFastaRecord(fasta, [sequence:true])
                .sequence == 'NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS\nNYITPVN\n'

        FastaSplitter.parseFastaRecord(fasta, [sequence:true, width: 20 ])
                .sequence == 'NLFVALYDFVASGDNTLSIT\nKGEKLRVLGYNHNGEWCEAQ\nTKNGQGWVPSNYITPVN\n'

        FastaSplitter.parseFastaRecord(fasta, [header:true])
                .header == 'gi|5524211|gb|AAD44166.1| cytochrome b'

        FastaSplitter.parseFastaRecord(fasta, [seqString:true])
                .seqString == 'NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVN'

        FastaSplitter.parseFastaRecord(fasta, [text:true])
                .text == fasta

        FastaSplitter.parseFastaRecord(fasta, [desc:true])
                .desc == 'cytochrome b'

    }



    def testSplitFasta () {

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
        def q = new FastaSplitter().options(each:{ count++; it }).target(fasta).channel()

        then:
        count == 2
        q.val == ">prot1\nLCLYTHIGRNIYYGS1\nEWIWGGFSVDKATLN\n"
        q.val == ">prot2\nLLILILLLLLLALLS\nGLMPFLHTSKHRSMM\nIENY\n"
        q.val == Channel.STOP

    }

    def testSplitFastaRecord() {

        given:
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

        when:
        def q = new FastaSplitter().options(record: [id:true, seqString:true]).target(fasta) .channel()

        then:
        q.val == [id:'1aboA', seqString: 'NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVN']
        q.val == [id:'1ycsB', seqString: 'KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGYVPRNLLGLYP']
        q.val == [id:'1pht', seqString: 'GYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIGWLNGYNETTGERGDFPGTYVEYIGRKKISP']
        q.val == Channel.STOP

    }

    def testSplitFastaFile () {

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
        def result = new FastaSplitter().options(by:2).target(fasta).list()

        then:
        result[0] == ">prot1\nAA\n>prot2\nBB\nCC\n"
        result[1] == ">prot3\nDD\n>prot4\nEE\nFF\nGG\n"
        result[2] == ">prot5\nLL\nNN\n"


        when:
        def result2 = new FastaSplitter()
                .options(record: [id: true, seqString: true], each:{ [ it.id, it.seqString.size() ]} )
                .target(fasta)
                .list()

        then:
        result2[0] == [ 'prot1', 2 ]
        result2[1] == [ 'prot2', 4 ]
        result2[2] == [ 'prot3', 2 ]
        result2[3] == [ 'prot4', 6 ]
        result2[4] == [ 'prot5', 4 ]
        result2.size() == 5

    }


    def testSplitWithLimit() {

        given:
        def fasta = '''
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
            .stripIndent().leftTrim()

        when:
        def result = new FastaSplitter(record:[id: true], limit: 3).target(fasta).list()
        then:
        result.size() == 3
        result[0] == [id: '1aboA']
        result[1] == [id: '1ycsB']
        result[2] == [id: '1pht']

    }


    def testSplitToFile() {

        given:
        def folder = TestHelper.createInMemTempDir()
        def fasta = '''
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
                .stripIndent().leftTrim()

        when:
        def result = new FastaSplitter().options(by: 2, file: folder).target(fasta).list()
        then:
        result.size() == 3
        result[0].text ==
            '''
            >1aboA
            NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS
            NYITPVN
            >1ycsB
            KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY
            VPRNLLGLYP
            '''
                    .stripIndent().leftTrim()


        result[1].text == '''
            >1pht
            GYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG
            WLNGYNETTGERGDFPGTYVEYIGRKKISP
            >1vie
            DRVRKKSGAAWQGQIVGWYCTNLTPEGYAVESEAHPGSVQIYPVAALERI
            N
            '''
                .stripIndent().leftTrim()

        result[2].text == '''
            >1ihvA
            NFRVYYRDSRDPVWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRD
            '''
                        .stripIndent().leftTrim()
    }

    def testSplitToFileByOne() {

        given:
        def folder = TestHelper.createInMemTempDir()
        def fasta = '''
            >1aboA
            NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS
            NYITPVN
            >1ycsB
            KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY
            VPRNLLGLYP
            >1pht
            GYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG
            WLNGYNETTGERGDFPGTYVEYIGRKKISP
            '''
                .stripIndent().leftTrim()

        when:
        def result = new FastaSplitter().options(file: folder).target(fasta).list()
        then:
        result.size() == 3
        result[0].text ==
            '''
            >1aboA
            NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS
            NYITPVN
            '''
                    .stripIndent().leftTrim()

        result[1].text == '''
            >1ycsB
            KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY
            VPRNLLGLYP
            '''
                .stripIndent().leftTrim()

        result[2].text == '''
            >1pht
            GYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG
            WLNGYNETTGERGDFPGTYVEYIGRKKISP
            '''
                        .stripIndent().leftTrim()
    }


    def testSplitRecordBy2() {

        given:
        def fasta = '''
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
                .stripIndent().leftTrim()

        when:
        def result = new FastaSplitter().options(record:[id: true], by:2).target(fasta).list()
        then:
        result.size() == 3
        result[0] == [[id: '1aboA'], [id: '1ycsB']]
        result[1] == [[id: '1pht'], [id: '1vie']]
        result[2] == [[id: '1ihvA']]

    }

    def 'should split by size' () {
        given:
        def fasta = '''
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
                .stripIndent().leftTrim()

        when:
        def result = new FastaSplitter().options(size: '150 B').target(fasta).list()

        then:
        result[0] == '''
                    >1aboA
                    NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS
                    NYITPVN
                    >1ycsB
                    KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY
                    VPRNLLGLYP
                    >1pht
                    GYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG
                    WLNGYNETTGERGDFPGTYVEYIGRKKISP
                    '''.stripIndent().leftTrim()

        result[1] == '''
                    >1vie
                    DRVRKKSGAAWQGQIVGWYCTNLTPEGYAVESEAHPGSVQIYPVAALERI
                    N
                    >1ihvA
                    NFRVYYRDSRDPVWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRD
                    '''.stripIndent().leftTrim()
    }

    def 'should fetch a record' () {

        given:
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

        when:
        def splitter = new FastaSplitter()
        def result = splitter.fetchRecord(new BufferedReader(new StringReader(fasta)))
        then:
        result == '''
                >prot1
                LCLYTHIGRNIYYGS1
                EWIWGGFSVDKATLN
                '''.stripIndent().leftTrim()

        splitter.counter.increment == 1
        splitter.counter.size == 1

        when:
        splitter = new FastaSplitter().options(size: '1MB')
        result = splitter.fetchRecord(new BufferedReader(new StringReader(fasta)))
        then:
        result == '''
                >prot1
                LCLYTHIGRNIYYGS1
                EWIWGGFSVDKATLN
                '''.stripIndent().leftTrim()

        splitter.counter.increment == 31
        splitter.counter.size == 1024 * 1024

    }

}
