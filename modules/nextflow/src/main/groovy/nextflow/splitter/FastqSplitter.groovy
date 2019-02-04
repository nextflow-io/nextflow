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

package nextflow.splitter

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.InheritConstructors
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.exception.StopSplitIterationException
/**
 * Split FASTQ formatted text content or files
 *
 * @link http://en.wikipedia.org/wiki/FASTQ_format
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Emilio Palumbo <emilio.palumbo@crg.eu>
 */
@Slf4j
@CompileStatic
@InheritConstructors
class FastqSplitter extends AbstractTextSplitter {

    /**
     * Wrapper class that holds the index of the element splitter
     *
     * See {@link nextflow.extension.SplitterMergeClosure}
     */
    @ToString
    @EqualsAndHashCode
    @TupleConstructor
    static class SplitIndex {
        int value
        String toString() { "SplitIndex($value)" }
    }

    private boolean processQualityField

    private boolean emitSplitIndex

    private volatile boolean emitted

    @Override
    protected Map<String,Object> validOptions() {
        def result = super.validOptions()
        result.record = [ Boolean, Map ]
        return result
    }

    FastqSplitter setEmitSplitIndex( boolean value ) {
        emitSplitIndex = value
        return this
    }

    @PackageScope
    def findSource( List tuple ) {
        def result = super.findSource(tuple)

        if( emitSplitIndex && into instanceof DataflowWriteChannel && !emitted ) {
            log.trace("Emitting split index: $elem")
            append(into,new SplitIndex(elem))
            emitted = true
        }

        return result
    }

    static Map recordToMap( String l1, String l2, String l3, String l4, Map fields ) {
        def result = [:]

        if( !fields || fields.containsKey('readHeader'))
            result.readHeader = l1.substring(1)

        if( !fields || fields.containsKey('readString'))
            result.readString = l2

        if( !fields || fields.containsKey('qualityHeader'))
            result.qualityHeader = l3.substring(1)

        if( !fields || fields.containsKey('qualityString'))
            result.qualityString = l4

        return result
    }


    def private StringBuilder buffer = new StringBuilder()

    private String errorMessage = "Invalid FASTQ format"


    String recordToText( String l1, String l2, String l3, String l4 ) {
        buffer.setLength(0)

        // read header
        buffer << l1 << '\n'
        // read string
        buffer << l2 << '\n'
        // quality header
        buffer << l3 << '\n'
        // quality string
        buffer << l4 << '\n'

        return buffer.toString()
    }


    @Override
    protected process( Reader targetObject )  {

        if( sourceFile )
            errorMessage += " for file: " + sourceFile

        super.process( targetObject )
    }

    @Override
    protected fetchRecord(BufferedReader reader) {

        def l1 = reader.readLine()
        def l2 = reader.readLine()
        def l3 = reader.readLine()
        def l4 = reader.readLine()

        if( !l1 || !l2 || !l3 || !l4 )
            return null

        if( !l1.startsWith('@') || !l3.startsWith('+') )
            throw new IllegalStateException(errorMessage)

        if( recordMode )
            return recordToMap(l1,l2,l3,l4, recordFields)

        if( processQualityField )
            return l4

        return recordToText(l1,l2,l3,l4)
    }

    /**
     * Retrieve the encoding quality score of the fastq file
     *
     * See http://en.wikipedia.org/wiki/FASTQ_format#Encoding
     *
     * @param quality A fastq quality string
     */
    def int qualityScore(Map opts=null) {

        if( opts ) options(opts)

        processQualityField = true

        int result = -1
        closure = { String quality ->
            result = detectQualityString(quality)
            if( result != -1 )
                throw new StopSplitIterationException()
        }
        apply()

        return result
    }


    static int detectQualityString( String quality ) {
        if( !quality )
            return -1

        for (int i=0; i<quality.size(); i++) {
            def c = (int)quality.charAt(i)
            if( c < 59 )
                return 33
            if( c > 74 )
                return 64
        }

        return -1
    }
}
