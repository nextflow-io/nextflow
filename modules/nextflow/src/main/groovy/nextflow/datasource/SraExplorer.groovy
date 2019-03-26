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

import java.nio.file.NoSuchFileException
import java.nio.file.Path

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.ToString
import groovy.util.logging.Slf4j
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
/**
 * Query NCBI SRA database and returns the retrieved FASTQs to the specified
 * target channel. Inspired to SRA-Explorer by Phil Ewels -- https://ewels.github.io/sra-explorer/
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class SraExplorer {

    static public Map PARAMS = [apiKey:String, cache: Boolean, max: Integer]

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
    private int entriesPerChunk = 1000
    private List<String> missing = new ArrayList<>()
    private Path cacheFolder

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

    protected Map getEnv() { System.getenv() }

    protected Path getCacheFolder() {
        if( cacheFolder )
            return cacheFolder
        cacheFolder = Const.APP_HOME_DIR.resolve('ncbi/sra')
        cacheFolder.createDirIfNotExists()
        return cacheFolder
    }

    protected String getConfigApiKey() {
        def session = Global.session as Session
        def result = session ?.config ?. navigate('ncbi.apiKey')
        if( !result )
            result = getEnv().get('NCBI_API_KEY')
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
        final text = new URL(url).getText()

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
        final text = new URL(url).getText()

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

    protected String readRunUrl(String acc) {
        final url = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?result=read_run&fields=fastq_ftp&accession=$acc"
        log.debug "SRA fetch ftp fastq url=$url"
        String result = new URL(url).text.trim()
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

        def lines = text.trim().readLines()
        def files = lines[1].split(';')
        def result = new ArrayList(files.size())
        for( def str : files ) {
            result.add( FileHelper.asPath("ftp://$str") )
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

}


