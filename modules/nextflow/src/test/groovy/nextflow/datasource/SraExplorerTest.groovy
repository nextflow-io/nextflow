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
        def f0 = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/004/SRR1448774/SRR1448774.fastq.gz' as Path
        def f1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908503/ERR908503_1.fastq.gz" as Path
        def f2= "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908503/ERR908503_2.fastq.gz" as Path

        when:
        def result = slurper.getFastqUrl('SRR1448774')
        then:
        1 * slurper.readRunFastqs('SRR1448774') >> RESP1
        result == f0

        when:
        result = slurper.getFastqUrl('ERR908503')
        then:
        1 * slurper.readRunFastqs('ERR908503') >> RESP2
        result == [f1, f2]
    }

    def 'should cache fastq_ftp' () {
        given:
        def CACHE_CONTENT = 'Hello'
        def folder = Files.createTempDirectory('test')
        def cache = folder.resolve('dir1/cache.txt')
        def slurper = Spy(SraExplorer)

        when:
        def result = slurper.readRunFastqs('ERR908503')
        then:
        1 * slurper.cachePath("ERR908503") >> cache
        1 * slurper.readRunUrl("ERR908503") >> CACHE_CONTENT
        result == CACHE_CONTENT
        cache.text == CACHE_CONTENT

        // should HIT the cache
        when:
        result = slurper.readRunFastqs('ERR908503')
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
}
