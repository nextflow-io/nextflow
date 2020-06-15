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

package nextflow.datasource

import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Specification

import java.nio.file.Files
import java.nio.file.Path

import nextflow.Const
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SraExplorerTest extends Specification {

    def 'should return query url' () {
        when:
        def slurper = new SraExplorer(apiKey: API)
        def result = slurper.getSearchUrl(TERM)
        then:
        result == EXPECTED

        where:
        TERM            | API   | EXPECTED
        'foo'           | null  | 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&usehistory=y&retmode=json&term=foo'
        'foo and bar'   | null  | 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&usehistory=y&retmode=json&term=foo+and+bar'
        'xxx'           | '1234'| 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&usehistory=y&retmode=json&api_key=1234&term=xxx'
    }

    def 'should return data url' () {
        given:
        def slurper = new SraExplorer(apiKey: API)
        expect:
        slurper.getFetchUrl(KEY,ENV,START,MAX) == EXPECTED

        where:
        KEY     | ENV    | START | MAX   | API   | EXPECTED
        '1'     | 'abc'  | 1     | 10    | null  | 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&retmode=json&query_key=1&WebEnv=abc&retstart=1&retmax=10'
        '2'     | 'xyz'  | 5     | 100   | '1234'| 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&retmode=json&query_key=2&WebEnv=xyz&retstart=5&retmax=100&api_key=1234'

    }

    def 'should parse runs xml' () {

        given:
        def slurper = new SraExplorer()
        def runs = "&lt;Run acc=\"SRR534293\" total_spots=\"148437796\" total_bases=\"29984434792\" load_done=\"true\" is_public=\"true\" cluster_name=\"public\" static_data_available=\"true\"/&gt;&lt;Run acc=\"SRR534294\" total_spots=\"141727117\" total_bases=\"28628877634\" load_done=\"true\" is_public=\"true\" cluster_name=\"public\" static_data_available=\"true\"/&gt;"

        when:
        def result = slurper.parseXml(runs)
        then:
        result.Run.size() == 2
        result.Run[0].@acc.toString() == 'SRR534293'
        result.Run[1].@acc.toString() == 'SRR534294'
        result.Run[0].attributes() == [acc:'SRR534293', total_spots:'148437796', total_bases:'29984434792', load_done:'true', is_public:'true', cluster_name:'public', static_data_available:'true']
        result.Run[1].attributes() == [acc:'SRR534294', total_spots:'141727117', total_bases:'28628877634', load_done:'true', is_public:'true', cluster_name:'public', static_data_available:'true']

    }


    def 'should parse exp xml' () {

/*
       <Summary>
          <Title>GSM981245: CSHL_RnaSeq_MCF-7_nucleus_longPolyA</Title>
          <Platform instrument_model="Illumina HiSeq 2000">ILLUMINA</Platform>
          <Statistics total_runs="2" total_spots="290164913" total_bases="58613312426" total_size="37016329841" load_done="true" cluster_name="public" />
       </Summary>
       <Submitter acc="SRA039973" center_name="GEO" contact_name="Gene Expression Omnibus (GEO), NCBI, NLM, NIH, htt" lab_name="" />
       <Experiment acc="SRX174314" ver="1" status="public" name="GSM981245: CSHL_RnaSeq_MCF-7_nucleus_longPolyA" />
       <Study acc="SRP007461" name="GSE30567: Long RNA-seq from ENCODE/Cold Spring Harbor Lab" />
       <Organism taxid="9606" ScientificName="Homo sapiens" />
       <Sample acc="SRS353539" name="" />
       <Instrument ILLUMINA="Illumina HiSeq 2000" />
       <Library_descriptor>
          <LIBRARY_NAME>GSM981245: CSHL_RnaSeq_MCF-7_nucleus_longPolyA</LIBRARY_NAME>
          <LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>
          <LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>
          <LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION>
          <LIBRARY_LAYOUT>
             <PAIRED />
          </LIBRARY_LAYOUT>
       </Library_descriptor>
       <Biosample>SAMN01103815</Biosample>
 */

        given:
        def slurper = new SraExplorer()
        def exp = "&lt;Summary&gt;&lt;Title&gt;GSM981245: CSHL_RnaSeq_MCF-7_nucleus_longPolyA&lt;/Title&gt;&lt;Platform instrument_model=\"Illumina HiSeq 2000\"&gt;ILLUMINA&lt;/Platform&gt;&lt;Statistics total_runs=\"2\" total_spots=\"290164913\" total_bases=\"58613312426\" total_size=\"37016329841\" load_done=\"true\" cluster_name=\"public\"/&gt;&lt;/Summary&gt;&lt;Submitter acc=\"SRA039973\" center_name=\"GEO\" contact_name=\"Gene Expression Omnibus (GEO), NCBI, NLM, NIH, htt\" lab_name=\"\"/&gt;&lt;Experiment acc=\"SRX174314\" ver=\"1\" status=\"public\" name=\"GSM981245: CSHL_RnaSeq_MCF-7_nucleus_longPolyA\"/&gt;&lt;Study acc=\"SRP007461\" name=\"GSE30567: Long RNA-seq from ENCODE/Cold Spring Harbor Lab\"/&gt;&lt;Organism taxid=\"9606\" ScientificName=\"Homo sapiens\"/&gt;&lt;Sample acc=\"SRS353539\" name=\"\"/&gt;&lt;Instrument ILLUMINA=\"Illumina HiSeq 2000\"/&gt;&lt;Library_descriptor&gt;&lt;LIBRARY_NAME&gt;GSM981245: CSHL_RnaSeq_MCF-7_nucleus_longPolyA&lt;/LIBRARY_NAME&gt;&lt;LIBRARY_STRATEGY&gt;RNA-Seq&lt;/LIBRARY_STRATEGY&gt;&lt;LIBRARY_SOURCE&gt;TRANSCRIPTOMIC&lt;/LIBRARY_SOURCE&gt;&lt;LIBRARY_SELECTION&gt;cDNA&lt;/LIBRARY_SELECTION&gt;&lt;LIBRARY_LAYOUT&gt;                 &lt;PAIRED/&gt;               &lt;/LIBRARY_LAYOUT&gt;&lt;/Library_descriptor&gt;&lt;Biosample&gt;SAMN01103815&lt;/Biosample&gt;  "

        when:
        def result = slurper.parseXml(exp)
        then:
        result.Summary.Title.text() == 'GSM981245: CSHL_RnaSeq_MCF-7_nucleus_longPolyA'

    }

    def 'should parse exp csv' () {

/*
       study_accession	secondary_study_accession	sample_accession	secondary_sample_accession	experiment_accession	run_accession	submission_accession	tax_id	scientific_name	instrument_platform	instrument_model	library_name	nominal_length	library_layout	library_strategy	library_source	library_selection	read_count	base_count	center_name	first_public	last_updated	experiment_title	study_title	study_alias	experiment_alias	run_alias	fastq_bytes	fastq_md5	fastq_ftp	fastq_aspera	fastq_galaxy	submitted_bytes	submitted_md5	submitted_ftp	submitted_aspera	submitted_galaxy	submitted_format	sra_bytes	sra_md5	sra_ftp	sra_aspera	sra_galaxy	cram_index_ftp	cram_index_aspera	cram_index_galaxy	sample_alias	broker_name	sample_title	nominal_sdev	first_created
       PRJNA30709	SRP007461	SAMN00634070	SRS214591	SRX082565	SRR307897	SRA039973	9606	Homo sapiens	ILLUMINA	Illumina Genome Analyzer II	GSM758559: CshlLong_RnaSeq_GM12878_cell_longPolyA		PAIRED	RNA-Seq	TRANSCRIPTOMIC	cDNA	32214376	4896585152	GEO	2011-07-28	2015-06-19	Illumina Genome Analyzer II sequencing; GSM758559: CshlLong_RnaSeq_GM12878_cell_longPolyA	GSE30567: Long RNA-seq from ENCODE/Cold Spring Harbor Lab	GSE30567	GSM758559: lab_RnaSeq_GM12878_cell_longPolyA	GSM758559_1	2392496788;2436856506	cdeb9d2450ca0bfc678968ac33f551a6;a1e9011d47d5b0462b9d50402a92dc6d	ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307897/SRR307897_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307897/SRR307897_2.fastq.gz	fasp.sra.ebi.ac.uk:/vol1/fastq/SRR307/SRR307897/SRR307897_1.fastq.gz;fasp.sra.ebi.ac.uk:/vol1/fastq/SRR307/SRR307897/SRR307897_2.fastq.gz	ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307897/SRR307897_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307897/SRR307897_2.fastq.gz							3499798854	2ee356e58e9765d6b26bde4b119b99be	ftp.sra.ebi.ac.uk/vol1/srr/SRR307/SRR307897	fasp.sra.ebi.ac.uk:/vol1/srr/SRR307/SRR307897	ftp.sra.ebi.ac.uk/vol1/srr/SRR307/SRR307897				GSM758559		CSHL_RnaSeq_GM12878_cell_longPolyA (superseded by GSE86658)		2011-07-28

 */

        given:
        def slurper = new SraExplorer()
        def exp = "study_accession\tsecondary_study_accession\tsample_accession\tsecondary_sample_accession\texperiment_accession\trun_accession\tsubmission_accession\ttax_id\tscientific_name\tinstrument_platform\tinstrument_model\tlibrary_name\tnominal_length\tlibrary_layout\tlibrary_strategy\tlibrary_source\tlibrary_selection\tread_count\tbase_count\tcenter_name\tfirst_public\tlast_updated\texperiment_title\tstudy_title\tstudy_alias\texperiment_alias\trun_alias\tfastq_bytes\tfastq_md5\tfastq_ftp\tfastq_aspera\tfastq_galaxy\tsubmitted_bytes\tsubmitted_md5\tsubmitted_ftp\tsubmitted_aspera\tsubmitted_galaxy\tsubmitted_format\tsra_bytes\tsra_md5\tsra_ftp\tsra_aspera\tsra_galaxy\tcram_index_ftp\tcram_index_aspera\tcram_index_galaxy\tsample_alias\tbroker_name\tsample_title\tnominal_sdev\tfirst_created\n" +
                "PRJNA30709\tSRP007461\tSAMN00634070\tSRS214591\tSRX082565\tSRR307897\tSRA039973\t9606\tHomo sapiens\tILLUMINA\tIllumina Genome Analyzer II\tGSM758559: CshlLong_RnaSeq_GM12878_cell_longPolyA\t\tPAIRED\tRNA-Seq\tTRANSCRIPTOMIC\tcDNA\t32214376\t4896585152\tGEO\t2011-07-28\t2015-06-19\tIllumina Genome Analyzer II sequencing; GSM758559: CshlLong_RnaSeq_GM12878_cell_longPolyA\tGSE30567: Long RNA-seq from ENCODE/Cold Spring Harbor Lab\tGSE30567\tGSM758559: lab_RnaSeq_GM12878_cell_longPolyA\tGSM758559_1\t2392496788;2436856506\tcdeb9d2450ca0bfc678968ac33f551a6;a1e9011d47d5b0462b9d50402a92dc6d\tftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307897/SRR307897_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307897/SRR307897_2.fastq.gz\tfasp.sra.ebi.ac.uk:/vol1/fastq/SRR307/SRR307897/SRR307897_1.fastq.gz;fasp.sra.ebi.ac.uk:/vol1/fastq/SRR307/SRR307897/SRR307897_2.fastq.gz\tftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307897/SRR307897_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307897/SRR307897_2.fastq.gz\t\t\t\t\t\t\t3499798854\t2ee356e58e9765d6b26bde4b119b99be\tftp.sra.ebi.ac.uk/vol1/srr/SRR307/SRR307897\tfasp.sra.ebi.ac.uk:/vol1/srr/SRR307/SRR307897\tftp.sra.ebi.ac.uk/vol1/srr/SRR307/SRR307897\t\t\t\tGSM758559\t\tCSHL_RnaSeq_GM12878_cell_longPolyA (superseded by GSE86658)\t\t2011-07-28\n"

        when:
        def result = slurper.parseCsv(exp)
        then:
        result['study_accession'] == 'PRJNA30709'
        result['first_created'] == '2011-07-28'
    }

    def 'should return ftp files for accession id' () {
        given:
        def RESP1 = '''
                fastq_ftp
                ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/004/SRR1448774/SRR1448774.fastq.gz
                '''.stripIndent()

        def RESP2 = '''
                fastq_ftp
                ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908503/ERR908503_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908503/ERR908503_2.fastq.gz
                '''.stripIndent()


        def slurper = Spy(SraExplorer)
        def SRAfields = null
        def f0 = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/004/SRR1448774/SRR1448774.fastq.gz' as Path
        def f1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908503/ERR908503_1.fastq.gz" as Path
        def f2= "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908503/ERR908503_2.fastq.gz" as Path

        when:
        def result = slurper.getFastqUrl('SRR1448774', SRAfields)
        then:
        1 * slurper.readRunFastqs('SRR1448774', SRAfields) >> RESP1
        result == f0

        when:
        result = slurper.getFastqUrl('ERR908503', SRAfields)
        then:
        1 * slurper.readRunFastqs('ERR908503', SRAfields) >> RESP2
        result == [f1, f2]
    }

    def 'should return ftp files for accession id plus queried fields' () {
        given:
        def slurper = new SraExplorer()
        def id = 'SRR1448774'

        expect:
        // slurper.getFastqUrl(id,[FIELD])[FIELD] == EXPECTED
        slurper.getFastqUrl(id,[FIELD]) == EXPECTED

        where:
        FIELD                                           | EXPECTED
        'study_accession'                               | [fastq_ftp:'ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/004/SRR1448774/SRR1448774.fastq.gz', study_accession:'PRJNA253315']
        'secondary_study_accession'                     | [fastq_ftp:'ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/004/SRR1448774/SRR1448774.fastq.gz', secondary_study_accession:'SRP043510']
        'first_created'                                 | [fastq_ftp:'ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/004/SRR1448774/SRR1448774.fastq.gz', first_created:'2015-06-05']
        'instrument_platform'                           | [fastq_ftp:'ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/004/SRR1448774/SRR1448774.fastq.gz', instrument_platform:'ILLUMINA']
        'instrument_model'                              | [fastq_ftp:'ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/004/SRR1448774/SRR1448774.fastq.gz', instrument_model:'Illumina HiSeq 2000']
        'submitted_md5'                                 | [fastq_ftp:'ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/004/SRR1448774/SRR1448774.fastq.gz', submitted_md5:null]
        'study_accession,submitted_md5,first_created'   | [fastq_ftp:'ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/004/SRR1448774/SRR1448774.fastq.gz', study_accession:'PRJNA253315', submitted_md5:null, first_created:'2015-06-05']
    }

    def 'should cache fastq_ftp' () {
        given:
        def CACHE_CONTENT = 'Hello'
        def SRAfields = null
        def folder = Files.createTempDirectory('test')
        def cache = folder.resolve('dir1/cache.txt')
        def slurper = Spy(SraExplorer)

        when:
        def result = slurper.readRunFastqs('ERR908503', SRAfields)
        then:
        1 * slurper.cachePath("ERR908503") >> cache
        1 * slurper.readRunUrl("ERR908503", SRAfields) >> CACHE_CONTENT
        result == CACHE_CONTENT
        cache.text == CACHE_CONTENT

        // should HIT the cache
        when:
        result = slurper.readRunFastqs('ERR908503', SRAfields)
        then:
        1 * slurper.cachePath("ERR908503") >> cache
        0 * slurper.readRunUrl(_) >> null
        result == CACHE_CONTENT

        cleanup:
        folder?.deleteDir()
    }

    def 'should return cache path' () {
        given:
        def slurper = new SraExplorer()
        expect:
        slurper.cachePath('SRP043510') == Const.APP_HOME_DIR.resolve('ncbi/sra/11/SRP043510.fastq_ftp.cache')
    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    @Requires({System.getenv('NCBI_API_KEY')})
    def 'should explore sra' () {
        given:
        def key = System.getenv('NCBI_API_KEY')
        def slurper = new SraExplorer(query: 'SRP043510', apiKey: key, maxResults: 10)
        when:
        def target = slurper.apply()
        then:
        target.count().val == 10
    }

    def 'should retrieve NCBI api env' () {
        given:
        def slurper = Spy(SraExplorer)
        when:
        def result = slurper.getConfigApiKey()
        then:
        1 * slurper.getEnv() >> [NCBI_API_KEY: '1bc']
        then:
        result == '1bc'
    }

    def 'should return SRA fields from a string' () {
        given:
        def slurper = new SraExplorer()

        expect:
        slurper.getFields(FIELDS) == EXPECTED

        where:
        FIELDS    | EXPECTED
        'a'       | ['a']
        'a,b'     | ['a','b']
        'a, b'    | ['a','b']
        ['a']     | ['a']
        ['a','b'] | ['a','b']


    }
}
