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

package nextflow.splitter
import java.nio.charset.Charset
import java.nio.charset.CharsetDecoder
import java.nio.file.Path

import com.google.common.hash.HashCode
import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Nextflow
import nextflow.exception.AbortOperationException
import nextflow.util.CacheHelper
import nextflow.util.CharsetHelper

import nextflow.splitter.TextFileCollector.CachePath

/**
 * Implements an abstract text splitting for the main types
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
@InheritConstructors
abstract class AbstractTextSplitter extends AbstractSplitter<Reader> {

    protected Charset charset = Charset.defaultCharset()

    Charset getCharset() { charset }

    protected boolean fileMode

    protected Path collectPath

    protected String collectName

    protected boolean compress

    private long itemsCount

    AbstractTextSplitter options(Map options) {
        super.options(options)

        // the charset used to parse the file
        charset = CharsetHelper.getCharset(options.charset)

        // collect chunks to temporary files
        if( options.file instanceof Boolean )
            fileMode = options.file as Boolean
        else if( options.file instanceof Path ) {
            collectPath = (Path)options.file
            fileMode = true
        }
        else if( options.file instanceof CharSequence ) {
            collectName = options.file.toString()
            fileMode = true
        }

        if( options.compress?.toString() == 'true' )
            compress = options.compress as boolean

        // file mode and record cannot be used at the same time
        if( fileMode && recordMode )
            throw new AbortOperationException("Parameters `file` and `record` conflict -- check operator `$operatorName`")

        if( !fileMode && compress )
            throw new AbortOperationException("Parameter `compress` requires also the use of parameter `file: true` -- check operator `$operatorName`")

        return this
    }

    /**
     * @return A map representing the valid options for the splitter. The map keys define the
     * accepted parameter names, the values the valid values for each of them.
     */
    @Override
    protected Map<String,Object> validOptions() {
        def result = super.validOptions()
        result.charset = [ Charset, Map, String ]
        result.file = [Boolean, Path, CharSequence]
        result.compress = [Boolean]
        return result
    }

    /**
     * Creates a {@link Reader} for the given object.
     *
     * @param obj An instance of {@link Reader}, {@link CharSequence}, {@link Path}, {@link File}
     *      {@link InputStream}, {@code char[]}
     * @return  A  {@link Reader} for the given object.
     * @throws IllegalArgumentException if the object specified is of a type not supported
     */
    protected Reader normalizeSource( obj ) {

        if( obj instanceof Reader )
            return (Reader)obj

        if( obj instanceof CharSequence )
            return new StringReader(obj.toString())

        if( obj instanceof Path )
            return newReader((Path)obj, charset)

        if( obj instanceof InputStream )
            return new InputStreamReader((InputStream)obj,charset)

        if( obj instanceof File )
            return newReader(obj.toPath(), charset)

        if( obj instanceof char[] )
            return new StringReader(new String((char[])obj))

        throw new IllegalArgumentException("Object of class '${obj.class.name}' does not support 'splitter' methods")

    }

    /**
     * Create a new reader object for the given path
     * @param path
     * @param charset
     * @return
     */
    protected Reader newReader( Path path, Charset charset ) {
        def source = newInputStream(path)
        CharsetDecoder decoder = charset.newDecoder()
        Reader reader = new InputStreamReader(source, decoder)
        new BufferedReader(reader)
    }

    /**
     * Wrap the target object with a {@link BufferedReader} filter stream
     * @param targetObject The target stream go be split
     * @return A {@link BufferedReader} instance
     */
    protected BufferedReader wrapReader( Reader targetObject ) {
        if( targetObject instanceof BufferedReader )
            return (BufferedReader)targetObject
        else
            return new BufferedReader(targetObject)
    }

    protected boolean isCollectorEnabled() {
        return (counter.isEnabled() || fileMode)
    }

    /**
     * Process the object to split
     *
     * @param targetObject
     * @param offset
     * @return
     */
    protected process( Reader targetObject ) {

        def result = null
        BufferedReader reader = wrapReader(targetObject)
        counter.reset() // <-- make sure to start
        itemsCount = 0
        try {
            while( true ) {
                // -- parse a record object
                final record = fetchRecord( reader )
                if( record == null )
                    break

                // -- apply the splitting logic for the fetched record
                result = processChunk( record )

                // -- check the limit of allowed rows has been reached
                if( limit>0 && ++itemsCount == limit )
                    break
            }

            // make sure to process collected entries
            if ( collector && collector.hasChunk() ) {
                result = invokeEachClosure(closure, collector.nextChunk())
            }
        }

        finally {
            reader.closeQuietly()
            if( collector instanceof Closeable )
                collector.closeQuietly()
        }

        return result
    }

    /**
     * Add the fetched record to the current collection and emit a new chunk
     *
     * @param record The read record object
     * @return A value returned by the user splitting closure
     */
    protected processChunk( record ) {

        def result = null

        if ( isCollectorEnabled() ) {
            // -- append to the list buffer
            collector.add(record)

            if( counter.isChunckComplete() ) {
                result = invokeEachClosure(closure, collector.nextChunk())
                counter.reset()
            }
        }
        else {
            result = invokeEachClosure(closure, record)
        }

        return result
    }

    /**
     * @return A {@link CollectorStrategy} object implementing a concrete
     * strategy according the user provided options
     */
    protected CollectorStrategy createCollector() {

        if( !isCollectorEnabled() )
            return null

        if( recordMode )
            return new ObjectListCollector()

        if( fileMode ) {
            def base = getCollectorBaseFile()
            return new TextFileCollector(base, charset, compress)
        }

        return new CharSequenceCollector()
    }

    protected String getCollectFileName() {
        if( collectName ) {
            return multiSplit ? "${collectName}_${elem}" : collectName
        }

        if( sourceFile ) {
            def fileName = sourceFile.getName()
            if( fileName.endsWith('.gz') )
                fileName = fileName[0..-4]

            return fileName
        }

        return multiSplit ? "chunk_$elem" : 'chunk'
    }

    /**
     * @return A path where file chunks are cached
     */
    @PackageScope
    CachePath getCollectorBaseFile () {

        final fileName = getCollectFileName()
        log.trace "Splitter collector file name: $fileName"

        CachePath result
        if( collectPath ) {
            final file = sourceFile ?: targetObj
            final path = collectPath.isDirectory() ? collectPath.resolve(fileName) : collectPath
            final keys = [Global?.session?.uniqueId, file, getCacheableOptions()]
            final hash = CacheHelper.hasher(keys).hash()
            result = new CachePath(path,hash)
        }

        else if( sourceFile ) {
            final path = Nextflow.cacheableFile( [sourceFile, getCacheableOptions()], fileName)
            result = new CachePath(path)
        }

        else {
            final path = Nextflow.cacheableFile( [targetObj, getCacheableOptions()], fileName )
            result = new CachePath(path)
        }

        log.debug "Splitter `$operatorName` collector path: $result"
        return result
    }

    /**
     * @return A map of the options that contribute to create unique key
     */
    protected Map getCacheableOptions() {
        def result = new HashMap( fOptionsMap )
        result.remove('into')
        result.remove('each')
        result.charset = charset.toString()
        return result
    }

    /**
     * Read and parse a new record from the reader object
     * @param reader A reader representing the object to split
     * @return A logical record read from underlying stream
     */
    abstract protected fetchRecord( BufferedReader reader )

}
