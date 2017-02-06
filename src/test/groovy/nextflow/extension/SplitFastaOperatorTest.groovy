/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
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


    def 'should split fasta in sequences'() {

        given:
        def sequences = Channel.from(fasta1).splitFasta()

        expect:
        with(sequences) {
            val == '>1aboA\nNLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPS\nNYITPVN\n'
            val == '>1ycsB\nKGVIYALWDYEPQNDDELPMKEGDCMTIIHREDEDEIEWWWARLNDKEGY\nVPRNLLGLYP\n'
            val == '>1pht\nGYQYRALYDYKKEREEDIDLHLGDILTVNKGSLVALGFSDGQEARPEEIG\nWLNGYNETTGERGDFPGTYVEYIGRKKISP\n'
            val == Channel.STOP
        }

    }


    def 'should split fasta in records' () {

        given:
        def records = Channel.from(fasta1, fasta2).splitFasta(record:[id:true])
        expect:
        records.val == [id:'1aboA']
        records.val == [id:'1ycsB']
        records.val == [id:'1pht']
        records.val == [id:'alpha123']
        records.val == Channel.STOP

    }

    def 'should split tuple in fasta records' () {

        given:
        def result = Channel.from( [fasta1, 'one'], [fasta2,'two'] ).splitFasta(record:[id:true]) { record, code ->
            [record.id, code]
        }

        expect:
        result.val == ['1aboA', 'one']
        result.val == ['1ycsB', 'one']
        result.val == ['1pht',  'one']
        result.val == ['alpha123', 'two']
        result.val == Channel.STOP

    }

    def 'should split fasta and forward result into the specified channel' () {

        given:
        def target = Channel.create()
        Channel.from(fasta1,fasta2).splitFasta(record:[id:true], into: target)

        expect:
        target.val == [id:'1aboA']
        target.val == [id:'1ycsB']
        target.val == [id:'1pht']
        target.val == [id:'alpha123']
        target.val == Channel.STOP
    }

}
