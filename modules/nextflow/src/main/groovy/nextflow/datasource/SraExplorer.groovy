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

package nextflow.datasource

import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedSupplier

import java.nio.file.NoSuchFileException
import java.nio.file.Path

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import groovy.xml.XmlParser
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.Const
import nextflow.Global
import nextflow.Session
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
import nextflow.util.CacheHelper
import nextflow.util.Duration

import java.time.temporal.ChronoUnit
import java.util.function.Predicate
import java.util.regex.Pattern

/**
 * Query NCBI SRA database and returns the retrieved FASTQs to the specified
 * target channel. Inspired to SRA-Explorer by Phil Ewels -- https://ewels.github.io/sra-explorer/
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class SraExplorer {

    static public Map PARAMS = [apiKey:[String,GString], cache: Boolean, max: Integer, protocol: ['ftp','http','https'], retryPolicy: Map]
    final static public List<Integer> RETRY_CODES = List.of(408, 429, 500, 502, 503, 504)
    final static private Pattern ERROR_PATTERN = ~/Server returned HTTP response code: (\d+) for URL.*/

    @ToString
    static class SearchRecord {
        int count
        int retmax
        int retstart
        String querykey
        String webenv
        String term
    }

    private Duration CACHE_MAX_TIME = Duration.of('30 days')

    private query
    private DataflowWriteChannel target
    private JsonSlurper jsonSlurper = new JsonSlurper()
    private XmlParser xmlParser = new XmlParser()
    private int emitCount
    private int maxResults = Integer.MAX_VALUE
    private int entriesPerChunk = 500
    private List<String> missing = new ArrayList<>()
    private Path cacheFolder
    private String protocol = 'ftp'
    private SraRetryConfig retryConfig = new SraRetryConfig()

    String apiKey
    boolean useCache = true

    SraExplorer() {

    }

    SraExplorer(DataflowWriteChannel target, Map opts) {
        this.target = target
        init(opts)
    }

    SraExplorer setQuery(query) {
        this.query = query
        return this
    }

    protected void init(Map opts) {
        if( opts.apiKey )
            apiKey = opts.apiKey
        if( opts.cache != null )
            useCache = opts.cache as boolean
        if( opts.max )
            maxResults = opts.max as int
        if( opts.protocol )
            protocol = opts.protocol as String
        if( opts.retryPolicy )
            retryConfig = new SraRetryConfig(opts.retryPolicy as Map)
    }

    DataflowWriteChannel apply() {
        if( target == null )
            target = new DataflowQueue()

        if( !apiKey )
            apiKey = getConfigApiKey()

        query0(query)

        if( missing )
            log.warn "Failed to retrieve fastq download URL for accessions: ${missing.join(',')}"

        target.bind(Channel.STOP)
        return target
    }

    protected Map env() {
        return System.getenv()
    }

    protected Map config() {
        final session = Global.session as Session
        return session.getConfig()
    }

    protected Path getCacheFolder() {
        if( cacheFolder )
            return cacheFolder
        cacheFolder = Const.APP_HOME_DIR.resolve('ncbi/sra')
        cacheFolder.createDirIfNotExists()
        return cacheFolder
    }

    protected String getConfigApiKey() {
        def result = config().navigate('ncbi.apiKey')
        if( !result )
            result = env().get('NCBI_API_KEY')
        if( !result )
            log.warn1("Define the NCBI_API_KEY env variable to use NCBI search service -- Read more https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/")
        return result
    }

    protected void query0( query ) {
        if( query instanceof List ) {
            for( def item in query ) {
                query0(item)
            }
            return
        }

        if( query instanceof String || query instanceof GString ) {
            query1(query.toString())
            return
        }

        throw new IllegalArgumentException("Not a valid query argument: $query [${query?.getClass()?.getName()}]")
    }

    protected void query1(String query) {

        def url = getSearchUrl(query)
        def result = makeSearch(url)
        int index = result.retstart ?: 0

        while( index < result.count && emitCount<maxResults ) {
            url = getFetchUrl(result.querykey, result.webenv, index, entriesPerChunk)
            def data = makeDataRequest(url)
            if( data.error )
                throw new IllegalArgumentException("Invalid request: $url -- Error: $data.error")
            parseDataResponse(data)
            index += entriesPerChunk
        }

    }

    protected String getSearchUrl(String term) {
        def url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&usehistory=y&retmode=json"
        if( apiKey )
            url += "&api_key=$apiKey"
        url += "&term=${URLEncoder.encode(term, "UTF-8")}"
        return url
    }

    protected Map makeDataRequest(String url) {
        log.debug "SRA data request url=$url"
        final text = runWithRetry(()->getTextFormUrl(url))

        log.trace "SRA data result:\n${pretty(text)?.indent()}"
        def response = jsonSlurper.parseText(text)

        if( response instanceof Map )
            return response

        throw new IllegalStateException("Invalid response object for URL $url")
    }

    private String pretty(String text) {
        try {
            return JsonOutput.prettyPrint(text)
        }
        catch( Exception e ) {
            return text
        }
    }

    protected void parseDataResponse( Map response ) {
        def ids = response.result.uids

        for( def key : ids ) {
            def runs = parseXml(response.result[key]?.runs)
            for( def run : runs ) {
                final acc = run['@acc']?.toString()
                final files = getFastqUrl(acc)
                if( acc && files ) {
                    target.bind( [acc, files] )
                }
                if( ++emitCount>= maxResults )
                    return
            }
        }
    }

    protected SearchRecord makeSearch(String url) {
        log.debug "SRA search url=$url"
        final text = runWithRetry(()-> getTextFormUrl(url))

        log.trace "SRA search result:\n${pretty(text)?.indent()}"
        final response = jsonSlurper.parseText(text)

        if( response instanceof Map && response.esearchresult instanceof Map ) {
            def search = (Map)response.esearchresult
            def result = new SearchRecord()
            result.count = search.count as Integer 
            result.retmax = search.retmax as Integer
            result.retstart = search.retstart as Integer
            result.querykey = search.querykey
            result.webenv = search.webenv
            result.term = search.querytranslation
            return result
        }

        throw new IllegalStateException("Invalid response object for URL $url")
    }

    protected Path cachePath( String acc ) {
        final bucket = CacheHelper.hasher(acc).hash().toString().substring(0,2)
        getCacheFolder().resolve("$bucket/${acc}.fastq_ftp.cache")
    }

    protected String readRunFastqs(String acc) {
        final Path cache = cachePath(acc)
        if( useCache ) try {
            def delta = System.currentTimeMillis() - FilesEx.lastModified(cache)
            if( delta < CACHE_MAX_TIME.millis )
                return cache.getText()
        }
        catch( NoSuchFileException e ) {
            // ok ignore it and download it from EBI
        }

        final result = readRunUrl(acc)
        if( useCache && result ) {
            // store in the cache
            cache.parent.createDirIfNotExists()
            cache.text = result
        }
        return result
    }

    protected static String getTextFormUrl(String url) {
        new URI(url).toURL().getText()
    }

    protected String readRunUrl(String acc) {
        final url = "https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run&fields=fastq_ftp&accession=$acc"
        log.debug "SRA fetch ftp fastq url=$url"
        String result = runWithRetry(() -> getTextFormUrl(url)).trim()
        log.trace "SRA fetch ftp fastq url result:\n${result?.indent()}"

        if( result.indexOf('\n')==-1 ) {
            log.debug "Invalid fastq ftp file format -- accession=$acc url=$url"
            missing << acc
            return null
        }

        return result
    }

    protected getFastqUrl(String acc) {
        def text = readRunFastqs(acc)
        if( !text )
            return

        // it returns tab separated like structures
        def lines = text.trim().readLines()
        // the first line contains the field names eg. "run_accession	fastq_ftp"
        def fields = lines[0].tokenize('\t')
        def index = fields.indexOf('fastq_ftp')
        if( index==-1 ) {
            log.warn "Unable to find SRA fastq_ftp field - Offending value: ${lines[0]}"
            return
        }
        // the second line holds one more tab separated values eg.
        // "SRR1448795	ftp.sra.ebi.ac.uk/.../SRR1448795_1.fastq.gz;ftp.sra.ebi.ac.uk/.../SRR1448795_2.fastq.gz"
        def value = lines[1].tokenize('\t')[index]
        def files = value.split(';')
        def result = new ArrayList(files.size())
        for( def str : files ) {
            result.add( FileHelper.asPath("$protocol://$str") )
        }

        return result.size()==1 ? result[0] : result
    }


    protected parseXml(String str) {
        if( !str )
            return null

        def xml = str.trim().replaceAll('&lt;','<').replaceAll('&gt;','>')
        xmlParser.parseText("<document>$xml</document>")
    }

    /**
     * NCBI search https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESummary
     * 
     * @param result
     * @return
     */
    protected String getFetchUrl(String key, String webenv, int retstart, int retmax) {

        def url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&retmode=json&query_key=${key}&WebEnv=${webenv}&retstart=$retstart&retmax=$retmax"
        if( apiKey )
            url += "&api_key=$apiKey"

        return url
    }

    /**
     * Creates a retry policy using the SRA retry configuration
     *
     * @param cond A predicate that determines when a retry should be triggered
     * @return The {@link dev.failsafe.RetryPolicy} instance
     */
    protected <T> RetryPolicy<T> retryPolicy(Predicate<? extends Throwable> cond) {
        final EventListener<ExecutionAttemptedEvent> listener = new EventListener<ExecutionAttemptedEvent>() {
            @Override
            void accept(ExecutionAttemptedEvent event) throws Throwable {
                log.debug("Retryable response error - attempt: ${event.attemptCount}; reason: ${event.lastFailure.message}")
            }
        }
        return RetryPolicy.<T>builder()
            .handleIf(cond)
            .withBackoff(retryConfig.delay.toMillis(), retryConfig.maxDelay.toMillis(), ChronoUnit.MILLIS)
            .withMaxAttempts(retryConfig.maxAttempts)
            .withJitter(retryConfig.jitter)
            .onRetry(listener)
            .build()
    }

    /**
     * Carry out the invocation of the specified action using a retry policy
     * when {@link java.io.IOException} is returned containing an error code.
     *
     * @param action A {@link dev.failsafe.function.CheckedSupplier} instance modeling the action to be performed in a safe manner
     * @return The result of the supplied action
     */
    protected <T> T runWithRetry(CheckedSupplier<T> action) {
        // define listener
        final listener = new EventListener<ExecutionAttemptedEvent>() {
            @Override
            void accept(ExecutionAttemptedEvent event) throws Throwable {
                log.debug("Retryable response error - attempt: ${event.attemptCount}; reason: ${event.lastFailure.message}")
            }
        }
        // define the retry condition
        final cond = new Predicate<? extends Throwable>() {
            @Override
            boolean test(Throwable t) {
                if( t instanceof IOException && containsErrorCodes(t.message, RETRY_CODES))
                    return true
                if(t.cause instanceof IOException && containsErrorCodes(t.cause.message, RETRY_CODES))
                    return true
                return false
            }
        }
        // create the retry policy
        def policy = retryPolicy(cond)
        // apply the action with
        return Failsafe.with(policy).get(action)
    }

    static boolean containsErrorCodes(String message, List<Integer> codes){
        def matcher = (message =~ ERROR_PATTERN)
        def httpCode = matcher ? matcher[0][1] as Integer : null
        return httpCode != null && codes.contains(httpCode)
    }

}


