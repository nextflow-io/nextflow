/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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

package nextflow.cloud.aws

import java.nio.file.Files
import java.nio.file.Paths
import java.nio.file.attribute.FileTime

import nextflow.util.Duration
import spock.lang.Specification
import spock.lang.Unroll
import test.TestHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
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
        def fields = parser.parseCsvLine('"alpha","beta",,,"delt\"a",')
        then:
        fields[0] == 'alpha'
        fields[1] == 'beta'
        fields[2] == null
        fields[3] == null
        fields[4] == 'delt\"a'
        fields[5] == null
        fields.size() == 6

        when:
        def line = '"6QB2UH376KYJD6MB","6QCMYABX3D","6QB2UH376KYJD6MB.6QCMYABX3D.6YS6EN2CT7","Reserved","USD 0.0 per Red Hat Enterprise Linux (Amazon VPC), t2.medium instance-hour (or partial hour)","2015-04-30","0","Inf","Hrs","0.0000000000","USD","1yr","All Upfront","standard","Compute Instance","AmazonEC2","Asia Pacific (Sydney)","AWS Region","t2.medium",,"General purpose","2","Intel Xeon Family","Up to 3.3 GHz","4 GiB","EBS only","Low to Moderate","32-bit or 64-bit",,,,,,,,"Shared",,"RHEL","No License required",,,,,,,,"APS2-BoxUsage:t2.medium","RunInstances:0010",,,,,,,,,,,,,,,"NA","Intel AVX; Intel Turbo",'
        def result = parser.parseCsvLine(line)
        then:
        result.size()==65
        result[0] == '6QB2UH376KYJD6MB'
        result[AmazonPriceReader.colPricePerUnit] == '0.0000000000'
        result[AmazonPriceReader.colCurrency] == 'USD'
        result[AmazonPriceReader.colProdFamily] == 'Compute Instance'
        result[AmazonPriceReader.colServiceCode] == 'AmazonEC2'
        result[AmazonPriceReader.colInstanceType] == 't2.medium'
        result[AmazonPriceReader.colCpu] == '2'
        result[AmazonPriceReader.colMem] == '4 GiB'
        result[AmazonPriceReader.colLocation]  == "Asia Pacific (Sydney)"

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
