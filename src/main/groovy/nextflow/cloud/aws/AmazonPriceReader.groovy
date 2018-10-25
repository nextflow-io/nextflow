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

    final static private int CACHE_MAX_DAYS = 7

    final static private int COL_CURRENCY = 10
    final static private int COL_PRODFAMILY = 14 // it must be `Compute Instance`
    final static private int COL_SERVICECODE = 15 // it must be `AmazonEC2`
    final static private int COL_LOCATION = 16 //
    final static private int COL_INSTANCETYPE = 18
    final static private int COL_CPU = 21
    final static private int COL_MEM = 24
    final static private int COL_STORAGE = 25
    final static private int COL_OS = 36

    final static private char DOUBLE_QUOTE = '"' as char
    final static private char COMMA = ',' as char

    final static private Map<String,String> REGIONS = new HashMap<>(15)

    static {
        REGIONS.'ap-south-1' = 'Asia Pacific (Mumbai)'
        REGIONS.'eu-west-1' = "EU (Ireland)"
        REGIONS.'eu-west-2' = "EU (London)"
        REGIONS.'eu-central-1' = "EU (Frankfurt)"
        REGIONS.'ap-southeast-1' = "Asia Pacific (Singapore)"
        REGIONS.'ap-southeast-2' = "Asia Pacific (Sydney)"
        REGIONS.'ap-northeast-1' = "Asia Pacific (Tokyo)"
        REGIONS.'ap-northeast-2' = "Asia Pacific (Seoul)"
        REGIONS.'us-east-1' = "US East (N. Virginia)"
        REGIONS.'us-west-1' = "US West (N. California)"
        REGIONS.'us-west-2' = "US West (Oregon)"
        REGIONS.'sa-east-1' = "South America (Sao Paulo)"
        REGIONS.'ca-central-1' = "Canada (Central)"

    }

    private Map<String,CloudInstanceType> TABLE

    private String region

    private int colProdFamily = COL_PRODFAMILY
    private int colServiceCode = COL_SERVICECODE
    private int colLocation = COL_LOCATION
    private int colInstanceType = COL_INSTANCETYPE
    private int colCpu = COL_CPU
    private int colMem = COL_MEM
    private int colStorage = COL_STORAGE
    private int colOS = COL_OS

    AmazonPriceReader( String region ) {
        this.region = region
        if( !REGIONS.keySet().contains(region) )
            throw new IllegalArgumentException("Unknown AWS region: `$region`")
    }

    CloudInstanceType getInstanceType( String typeId ) {
        return getInstanceTypeTable().get(typeId)
    }

    @Lazy
    private Path CACHE_FOLDER = {
        def result = APP_HOME_DIR.resolve('aws')
        result.mkdirs()
        return result
    }()

    private synchronized Map<String,CloudInstanceType> getInstanceTypeTable() {
        if( TABLE == null ) {
            TABLE = getCachedTable()
        }
        return TABLE
    }

    private Map<String,CloudInstanceType> getCachedTable() {

        final path = CACHE_FOLDER.resolve("ec2-instance-types-${region}.bin")
        if( existsAndNotExpired(path)) {
            return (Map<String,CloudInstanceType>)KryoHelper.deserialize(path)
        }

        final result = parse(ENDPOINT)
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
        def priceFile = CACHE_FOLDER.resolve('ec2-pricing.csv')
        if( existsAndNotExpired(priceFile) ) {
            return priceFile
        }

        log.info "Fetching EC2 prices (it can take a few seconds depending your internet connection) .."
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

        int error = 0
        final location = REGIONS[ region ]
        final map = new HashMap<>()
        // data starts at 7th row
        final reader = cachedUrlBufferedReader(endpoint)
        try {
            int c=0
            String line
            while( (line=reader.readLine()) != null ) {
                try {
                    c++
                    if (c < 7) {
                        if (c == 6) parseCsvHeader(line)
                        continue
                    }

                    def tkns = parseCsvLine(line)
                    if (tkns.size() < 64) continue
                    if (tkns[colProdFamily] != 'Compute Instance') continue
                    if (tkns[colServiceCode] != "AmazonEC2") continue
                    if (tkns[colLocation] != location) continue
                    if (tkns[colOS] != 'Linux') continue
                    List storage = parseStorage(tkns[colStorage])

                    final instanceType = tkns[colInstanceType]
                    if (map.get(instanceType))
                        continue

                    if (!tkns[colCpu]) {
                        log.debug "Missing cpu number for instance type: $instanceType -- offending entry: $line"
                        error++
                        continue
                    }

                    if (!tkns[colMem]) {
                        log.debug "Missing memory value for instance type: $instanceType -- offending entry: $line"
                        error++
                        continue
                    }

                    def entry = new CloudInstanceType(
                            id: instanceType,
                            cpus: tkns[colCpu] as int,
                            memory: parseMem(tkns[colMem]),
                            disk: storage[0],
                            numOfDisks: storage[1]
                    )

                    map.put(entry.id, entry)
                    log.debug "map store >> entry=$entry"
                }
                catch( Exception e ) {
                    log.debug "Unexptected AWS price entry -- offending line: $line"
                    error++
                }
            }
        }
        finally {
            reader.close()
        }

        if( error ) {
            log.warn "One or more errors have been detected parsing EC2 prices -- See log file for details"
        }

        return map
    }

    private void parseCsvHeader( String line ) {
        final header = parseCsvLine(line)

        colProdFamily   = ndx(header, "Product Family", COL_PRODFAMILY)
        colServiceCode  = ndx(header, "serviceCode", COL_SERVICECODE)
        colLocation     = ndx(header, "Location", COL_LOCATION)
        colInstanceType = ndx(header,"Instance Type", COL_INSTANCETYPE)
        colCpu          = ndx(header,"vCPU", COL_CPU)
        colMem          = ndx(header,"Memory", COL_MEM)
        colStorage      = ndx(header,"Storage", COL_STORAGE)
        colOS           = ndx(header,"Operating System",COL_OS)

        log.debug "AWS csv fields format: ProductFamily=$colProdFamily; ServiceCode=$colServiceCode; Location=$colLocation; InstanceType=$colInstanceType; Cpu=$colCpu; Mem=$colMem; Storage=$colStorage; OS=$colOS"
    }

    private int ndx(List<String> header, String field, int defValue) {
        def p = header.indexOf(field)
        if( p==-1 ) {
            log.warn "Missing field `$field` in AWS csv price file"
            return defValue
        }
        return p
    }

    @PackageScope
    List<String> parseCsvLine( String line ) {
        def result = []
        while( line != null ) {
            if( !line ) {
                result.add(null)
                break
            }
            else if( line.startsWith('"') ) {
                line = readQuotedValue(line, result)
            }
            else {
                line = readSimpleValue(line, result)
            }
        }
        return result
    }

    private String readSimpleValue(String line, List<String> result) {
        def p = line.indexOf( (int)COMMA )
        if( p == -1 ) {
            result.add(line)
            return null
        }
        else {
            result.add(line.substring(0,p) ?: null)
            return line.substring(p+1)
        }
    }

    private String readQuotedValue(String line, List<String> result) {
        def strip = line.substring(1)
        def p = strip.indexOf( (int)DOUBLE_QUOTE )
        if( p == -1 )
            throw new IllegalStateException("Missing double-quote termination in CSV value -- offending line: $line")
        result.add( strip.substring(0,p) )

        def next = p+1
        if( next<strip.size() ) {
            if( strip.charAt(next) != COMMA )
                throw new IllegalStateException("Invalid CSV value -- offending line: $line")
            return strip.substring(next+1)
        }
        else {
            return null
        }
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
