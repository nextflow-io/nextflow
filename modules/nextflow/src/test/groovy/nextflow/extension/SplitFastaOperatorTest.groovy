/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.extension

import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Channel
import nextflow.Session
import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Timeout

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class SplitFastaOperatorTest extends Specification {

    def setupSpec() {
        new Session()
    }

    @Shared
    def fasta1 = """\
                >1aboA
                NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS
                NYITPVN
                >1ycsB
                KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY
                VPRNLLGLYP
                ; comment
                >1pht
                GYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG
                WLNGYNETTGERGDFPGTYVEYIGRKKISP
                """.stripIndent()

    @Shared
    def fasta2 = '''
            >alpha123
            WLNGYNETTGERGDFPGTYVEYIGRKKISP
            VPRNLLGLYP
            '''

    @Shared
    def fasta3 = """\
                >1pzqA
                GYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG
                WLNGYNETTGERGDFPGTYVEYIGRKKISP
                >1xdtB
                KGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY
                VPRNLLGLYP
                >1bcdB
                NLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS
                NYITPVN
                """.stripIndent()

    def 'should split fasta in sequences'() {
        given:
        def sequences = Channel.of(fasta1).splitFasta() as DataflowReadChannel

        expect:
        sequences.unwrap() == '>1aboA\nNLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS\nNYITPVN\n'
        sequences.unwrap() == '>1ycsB\nKGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY\nVPRNLLGLYP\n'
        sequences.unwrap() == '>1pht\nGYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG\nWLNGYNETTGERGDFPGTYVEYIGRKKISP\n'
        sequences.unwrap() == Channel.STOP
    }

    def 'should split fasta in records' () {
        given:
        def records = Channel.of(fasta1, fasta2).splitFasta(record:[id:true])
        expect:
        records.unwrap() == [id:'1aboA']
        records.unwrap() == [id:'1ycsB']
        records.unwrap() == [id:'1pht']
        records.unwrap() == [id:'alpha123']
        records.unwrap() == Channel.STOP
    }

    def 'should split tuple in fasta records' () {
        given:
        def result = Channel
                .of( [fasta1, 'one'], [fasta2,'two'] )
                .splitFasta(record:[id:true])
                .map{ record, code -> [record.id, code] }

        expect:
        result.unwrap() == ['1aboA', 'one']
        result.unwrap() == ['1ycsB', 'one']
        result.unwrap() == ['1pht',  'one']
        result.unwrap() == ['alpha123', 'two']
        result.unwrap() == Channel.STOP
    }

    def 'should split fasta and forward result into the specified channel' () {
        given:
        def target = CH.queue()
        Channel.of(fasta1,fasta2).splitFasta(record:[id:true], into: target)

        expect:
        target.unwrap() == [id:'1aboA']
        target.unwrap() == [id:'1ycsB']
        target.unwrap() == [id:'1pht']
        target.unwrap() == [id:'alpha123']
        target.unwrap() == Channel.STOP
    }

    def 'should apply count on multiple entries'() {
        given:
        def F1 = '''
            >1
            AAA
            >2
            BBB
            >3
            CCC
            '''
                .stripIndent().trim()
        def F3 = '''
            >1
            EEE
            >2
            FFF
            >3
            GGG
            '''
                .stripIndent().trim()

        def target = CH.queue()

        when:
        Channel.of(F1,F3).splitFasta(by:2, into: target)
        then:
        target.unwrap() == '>1\nAAA\n>2\nBBB\n'
        target.unwrap() == '>3\nCCC\n'
        target.unwrap() == '>1\nEEE\n>2\nFFF\n'
        target.unwrap() == '>3\nGGG\n'
    }

    def 'should apply count on multiple entries with a limit'() {
        given:
        def F1 = '''
            >1
            AAA
            >2
            BBB
            >3
            CCC
            >4
            DDD
            >5
            XXX
            '''
                .stripIndent().trim()
        def F3 = '''
            >1
            EEE
            >2
            FFF
            >3
            GGG
            >4
            HHH
            >5
            YYY
            '''
                .stripIndent().trim()

        def target = CH.queue()

        when:
        Channel.of(F1,F3).splitFasta(by:2, limit:4, into: target)
        then:
        target.unwrap() == '>1\nAAA\n>2\nBBB\n'
        target.unwrap() == '>3\nCCC\n>4\nDDD\n'
        target.unwrap() == '>1\nEEE\n>2\nFFF\n'
        target.unwrap() == '>3\nGGG\n>4\nHHH\n'
    }
}
