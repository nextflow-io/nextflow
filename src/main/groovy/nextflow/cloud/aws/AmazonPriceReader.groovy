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
import static nextflow.Const.APP_HOME_DIR

import java.nio.charset.Charset
import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.cloud.types.CloudInstanceType
import nextflow.util.Duration
import nextflow.util.KryoHelper
import nextflow.util.MemoryUnit
/**
 * Download and parse Amazon AWS cloud price CSV file
 *
 * See  http://docs.aws.amazon.com/awsaccountbilling/latest/aboutv2/price-changes.html
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AmazonPriceReader {

    final static public String ENDPOINT = 'https://pricing.us-east-1.amazonaws.com/offers/v1.0/aws/AmazonEC2/current/index.csv'

    final static public int CACHE_MAX_DAYS = 7

    final static private int colPricePerUnit = 9
    final static private int colCurrency = 10
    final static private int colProdFamily = 13 // it must be `Compute Instance`
    final static private int colServiceCode = 14 // it must be `AmazonEC2`
    final static private int colLocation = 15 //
    final static private int colInstanceType = 17
    final static private int colCurrentGen = 18
    final static private int colCpu = 20
    final static private int colMem = 23
    final static private int colStorage = 24
    final static private int colOS = 36

    final static private char DOUBLE_QUOTE = '"' as char
    final static private char COMMA = ',' as char

    final static private Map<String,String> REGIONS = new HashMap<>(15)

    static {
        REGIONS.'ap-south-1' = 'Asia Pacific (Mumbai)'
        REGIONS.'eu-west-1' = "EU (Ireland)"
        REGIONS.'eu-central-1' = "EU (Frankfurt)"
        REGIONS.'ap-southeast-1' = "Asia Pacific (Singapore)"
        REGIONS.'ap-southeast-2' = "Asia Pacific (Sydney)"
        REGIONS.'ap-northeast-1' = "Asia Pacific (Tokyo)"
        REGIONS.'ap-northeast-2' = "Asia Pacific (Seoul)"
        REGIONS.'us-east-1' = "US East (N. Virginia)"
        REGIONS.'us-west-1' = "US West (N. California)"
        REGIONS.'us-west-2' = "US West (Oregon)"
        REGIONS.'sa-east-1' = "South America (Sao Paulo)"
    }

    private Map<String,CloudInstanceType> TYPES_CACHE

    private final Path EC2_TYPES_FILE

    private final Path EC2_PRICING_FILE

    private String region

    AmazonPriceReader( String region ) {
        this.region = region
        if( !REGIONS.keySet().contains(region) )
            throw new IllegalArgumentException("Unknown AWS region: `$region`")
        EC2_PRICING_FILE = CACHE_FOLDER.resolve('ec2-pricing.csv')
        EC2_TYPES_FILE = CACHE_FOLDER.resolve("ec2-instance-types-${region}.bin")
    }

    @Lazy
    static public Path CACHE_FOLDER = {
        def result = APP_HOME_DIR.resolve('aws')
        result.mkdirs()
        return result
    }()

    synchronized Map<String,CloudInstanceType> getInstanceTypeTable() {
        if( TYPES_CACHE == null ) {
            TYPES_CACHE = getCachedTable()
        }
        return TYPES_CACHE
    }

    private Map<String,CloudInstanceType> getCachedTable() {

        def path = EC2_TYPES_FILE
        if( existsAndNotExpired(path)) {
            return (Map<String,CloudInstanceType>)KryoHelper.deserialize(path)
        }

        def result = parse(ENDPOINT)
        KryoHelper.serialize(result, path)
        return result
    }

    private boolean existsAndNotExpired( Path path ) {
        if( !path.exists() ) return false
        !evict(path, Duration.of("$CACHE_MAX_DAYS day"))
    }

    private boolean evict( Path path, Duration maxAge, long now = System.currentTimeMillis() ) {
        if( path.exists() ) {
            final age = Duration.of(now-path.lastModified())
            if( age > maxAge ) {
                path.delete()
                return true
            }
        }
        return false
    }

    synchronized private Path cachedUrl(String url) {
        def priceFile = EC2_PRICING_FILE
        if( existsAndNotExpired(priceFile) ) {
            return priceFile
        }

        // download the file
        def _in = new URL(url).openStream()
        try {
            Files.copy(_in, priceFile)
        }
        finally {
            _in.closeQuietly()
        }

        return priceFile
    }

    private BufferedReader cachedUrlBufferedReader(String url) {
        Files.newBufferedReader( cachedUrl(url), Charset.forName('utf-8') )
    }

    /**
     * Download and parse AWS price file
     *
     * @param endpoint The url from where download the CSV price file
     * @return A map of the extracted
     */
    Map<String,CloudInstanceType> parse(String endpoint) {

        final location = REGIONS[ region ]
        final map = new HashMap<>()
        // data starts at 7th row
        final reader = cachedUrlBufferedReader(endpoint)
        try {
            int c=0
            String line
            while( (line=reader.readLine()) != null ) {
                if( c++ < 7 ) continue
                def tkns = parseCsvLine(line)
                if( tkns.size() < 64 ) continue
                if( tkns[colProdFamily] != 'Compute Instance') continue
                if( tkns[colServiceCode] != "AmazonEC2" ) continue
                if( tkns[colLocation] != location) continue
                if( tkns[colOS] != 'Linux' ) continue
                //if( tkns[colCurrentGen] != 'Yes' ) continue

                List storage = parseStorage(tkns[colStorage])

                if( map.get(tkns[colInstanceType]) )
                    continue

                def entry = new CloudInstanceType(
                        id: tkns[colInstanceType],
                        cpus: tkns[colCpu] as int ,
                        memory: parseMem(tkns[colMem]),
                        disk: storage[0],
                        numOfDisks: storage[1]
                )

                map.put(entry.id, entry)
            }
        }
        finally {
            reader.close()
        }

        return map
    }

    @PackageScope
    List<String> parseCsvLine( String line ) {
        def result = []
        boolean open = false
        int start = -1
        int i=0
        while( i<line.size() ) {

            def ch = line.charAt(i)
            if( ch == DOUBLE_QUOTE ) {
                if( !open ) {
                    // open a new field
                    open = true
                    start = i
                }
                else if( i+1 == line.size() || line.charAt(i+1) == COMMA ) {
                    // close the field
                    open = false
                    result << line.substring(start+1, i)
                }
            }
            else if( ch == COMMA ) {
                if( i+1 == line.size() || line.charAt(i+1) == COMMA ) {
                    // empty field
                    result << null
                }
            }

            i++
        }
        return result
    }

    private MemoryUnit parseMem( String str ) {
        new MemoryUnit(str.replace(',','').replace('i','').toUpperCase())
    }

    private static List EMPTY = [ MemoryUnit.ZERO , 0]

    private List parseStorage( String str ) {
        if( !str ) return EMPTY

        def item = str.tokenize(' ')
        if( item[0].isNumber() ) {
            if(item[1] == 'x' ) {
                return [ parseMem(item[2] +' GB'), item[0] as int ]
            }
            else {
                return [ parseMem(str), 1 ]
            }
        }
        else {
            return EMPTY
        }
    }


}
