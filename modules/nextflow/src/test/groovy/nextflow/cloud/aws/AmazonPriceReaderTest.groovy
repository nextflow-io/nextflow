/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.cloud.aws

import spock.lang.IgnoreIf
import spock.lang.Specification
import spock.lang.Unroll

import java.nio.file.Files
import java.nio.file.Paths
import java.nio.file.attribute.FileTime

import nextflow.util.Duration
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@IgnoreIf({System.getenv('NXF_QUICK')})
class AmazonPriceReaderTest extends Specification {

    static String REGION = 'eu-west-1'

    def setupSpec() {
        // delete the the pricing cached file
        def home = System.getProperty('user.home')
        def cache = Paths.get("$home/.nextflow/aws")
        cache.deleteDir()
    }


    @Unroll
    def 'parse price api #e_type' () {

        given:
        def result = new AmazonPriceReader(REGION).parse(AmazonPriceReader.ENDPOINT)

        expect:
        with(result.get(e_type)) {
            id == e_type
            cpus == e_cpus
            memory.toString() == e_mem
            disk.toString() == e_disk
            numOfDisks == e_num
        }

        where:
        e_type      | e_cpus | e_mem    | e_disk    | e_num
        't2.micro'  | 1      | '1 GB'   | '0'       | 0
        't2.medium' | 2      | '4 GB'   | '0'       | 0
        'm4.xlarge' | 4      | '16 GB'  | '0'       | 0
        'm3.xlarge' | 4      | '15 GB'  | '40 GB'   | 2
        'c4.2xlarge'| 8      | '15 GB'  | '0'       | 0
        'c3.8xlarge'| 32     | '60 GB'  | '320 GB'  | 2

    }



    def 'should parse a csv line' () {

        given:
        def parser = new AmazonPriceReader('eu-west-1')

        when:
        def fields = parser.parseCsvLine('"alpha","beta",,,"delta",')
        then:
        fields[0] == 'alpha'
        fields[1] == 'beta'
        fields[2] == null
        fields[3] == null
        fields[4] == 'delta'
        fields[5] == null
        fields.size() == 6

        when:
        fields = parser.parseCsvLine('"alpha",beta,"delta"')
        then:
        fields[0] == 'alpha'
        fields[1] == 'beta'
        fields[2] == 'delta'
        fields.size() == 3

        when:
        def line = '"6QB2UH376KYJD6MB","6QCMYABX3D","6QB2UH376KYJD6MB.6QCMYABX3D.6YS6EN2CT7","Reserved","USD 0.0 per Red Hat Enterprise Linux (Amazon VPC), t2.medium instance-hour (or partial hour)","2015-04-30","0","Inf","Hrs","0.0000000000","USD","1yr","All Upfront","standard","Compute Instance","AmazonEC2","Asia Pacific (Sydney)","AWS Region","t2.medium",,"General purpose","2","Intel Xeon Family","Up to 3.3 GHz","4 GiB","EBS only","Low to Moderate","32-bit or 64-bit",,,,,,,,"Shared",,"RHEL","No License required",,,,,,,,"APS2-BoxUsage:t2.medium","RunInstances:0010",,,,,,,,,,,,,,,"NA","Intel AVX; Intel Turbo",'
        def result = parser.parseCsvLine(line)
        then:
        result.size()==65
        result[0] == '6QB2UH376KYJD6MB'
        result[AmazonPriceReader.COL_CURRENCY] == 'USD'
        result[AmazonPriceReader.COL_PRODFAMILY] == 'Compute Instance'
        result[AmazonPriceReader.COL_SERVICECODE] == 'AmazonEC2'
        result[AmazonPriceReader.COL_INSTANCETYPE] == 't2.medium'
        result[AmazonPriceReader.COL_CPU] == '2'
        result[AmazonPriceReader.COL_MEM] == '4 GiB'
        result[AmazonPriceReader.COL_LOCATION]  == "Asia Pacific (Sydney)"

        when:
        line = '6S5B59JFU3B2M3XM,4NA7Y494T4,6S5B59JFU3B2M3XM.4NA7Y494T4.6YS6EN2CT7,Reserved,"Red Hat Enterprise Linux (Amazon VPC), x1.16xlarge reserved instance applied",2017-04-30,0,Inf,Hrs,5.0620000000,USD,1yr,No Upfront,standard,Compute Instance,AmazonEC2,EU (Ireland),AWS Region,x1.16xlarge,Yes,Memory optimized,64,High Frequency Intel Xeon E7-8880 v3 (Haswell),2.3 GHz,976 Gib,"1 x 1,920",High,64-bit,,,,,,,,Shared,RHEL,No License required,,,,,,,,EU-BoxUsage:x1.16xlarge,RunInstances:0010,,,124.5,Yes,,,,,,,,,,,Yes,Yes,Yes,128,,NA,,Amazon Elastic Compute Cloud'
        result = parser.parseCsvLine(line)
        then:
        result[0] == '6S5B59JFU3B2M3XM'
        result[AmazonPriceReader.COL_CURRENCY] == 'USD'
        result[AmazonPriceReader.COL_PRODFAMILY] == 'Compute Instance'
        result[AmazonPriceReader.COL_SERVICECODE] == 'AmazonEC2'
        result[AmazonPriceReader.COL_INSTANCETYPE] == 'x1.16xlarge'
        result[AmazonPriceReader.COL_CPU] == '64'
        result[AmazonPriceReader.COL_MEM] == '976 Gib'

        when:
        line = '"QDG74QMT7ZXABFC9","JRTCKXETXF","QDG74QMT7ZXABFC9.JRTCKXETXF.6YS6EN2CT7","OnDemand","$0.000 per Linux m4.large Dedicated Host Instance hour","2017-08-01","0","Inf","Hrs","0.0000000000","USD",,,,"Compute Instance","AmazonEC2","EU (Ireland)","AWS Region","m4.large","Yes","General purpose","2","Intel Xeon E5-2676 v3 (Haswell)","2.4 GHz","8 GiB","EBS only","Moderate","64-bit",,,,,,,,"Host","Linux","No License required",,,,,,,,"EU-HostBoxUsage:m4.large","RunInstances",,"450 Mbps","6.5","Yes",,,,,,,,,,,,,,"4",,"NA","Intel AVX; Intel AVX2; Intel Turbo","Amazon Elastic Compute Cloud"'
        result = parser.parseCsvLine(line)
        then:
        result[0] == 'QDG74QMT7ZXABFC9'
        result[AmazonPriceReader.COL_CURRENCY] == 'USD'
        result[AmazonPriceReader.COL_PRODFAMILY] == 'Compute Instance'
        result[AmazonPriceReader.COL_SERVICECODE] == 'AmazonEC2'
        result[AmazonPriceReader.COL_INSTANCETYPE] == 'm4.large'
        result[AmazonPriceReader.COL_CPU] == '2'
        result[AmazonPriceReader.COL_MEM] == '8 GiB'
        result[AmazonPriceReader.COL_OS] == 'Linux'
        result[AmazonPriceReader.COL_LOCATION]  == "EU (Ireland)"


    }

    def 'should validate path evict' ( ) {
        given:
        def parser = new AmazonPriceReader('eu-west-1')
        def now = System.currentTimeMillis()
        def _30secs = Duration.of('30 sec')

        def path1 = Files.createFile(TestHelper.createInMemTempFile('test'))
        Files.setLastModifiedTime(path1, FileTime.fromMillis(now - 60_000));

        def path2 = Files.createFile(TestHelper.createInMemTempFile('test'))
        Files.setLastModifiedTime(path2, FileTime.fromMillis(now - 10_000));


        expect:
        parser.evict(path1, _30secs, now)
        !path1.exists()

        !parser.evict(path2, _30secs, now)
        path2.exists()
    }


    def 'should validate exit and not expire method' () {
        given:
        def parser = new AmazonPriceReader('eu-west-1')
        def now = System.currentTimeMillis()

        def path1 = Files.createFile(TestHelper.createInMemTempFile('test'))
        Files.setLastModifiedTime(path1, FileTime.fromMillis(now));

        def path2 = Files.createFile(TestHelper.createInMemTempFile('test'))
        Files.setLastModifiedTime(path2, FileTime.fromMillis(now - AmazonPriceReader.CACHE_MAX_DAYS * 24 * 60 * 60 * 1000 * 2));

        expect:
        parser.existsAndNotExpired(path1)
        path1.exists()

        !parser.existsAndNotExpired(path2)
        !path2.exists()

        !parser.existsAndNotExpired(Paths.get('unknown'))
        !path2.exists()

    }
}
