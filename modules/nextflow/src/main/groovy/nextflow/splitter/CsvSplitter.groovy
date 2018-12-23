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
import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import org.apache.commons.lang.StringUtils
/**
 * Split a CSV file in records
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
@InheritConstructors
class CsvSplitter extends AbstractTextSplitter {

    /**
     * The CVS separator
     */
    protected String sep = ','

    /**
     * The character used to wrap values containing a separator char
     */
    protected String quote

    /**
     * When {@code true} remove leading and trailing blanks from a value
     */
    protected boolean stripBlanks

    /**
     * When {@code true} parse the first line as the columns names
     */
    protected boolean firstLineAsHeader

    /**
     * Defines the columns names of the resulting records
     */
    protected List<String> columnsHeader

    /**
     * The number of lines to skip
     */
    protected int skipLines = 0

    /**
     * Set the splitter options by specifying a map of named parameters.
     * Valid parameters are:
     * <li>{@code sep}
     * <li>{@code strip}
     * <li>{@code header}
     * <li>{@code quote}
     * <li>{@code skip}
     *
     * @param options
     * @return The splitter instance itself
     */
    CsvSplitter options(Map options) {

        super.options(options)

        // the separator character
        if( options.sep )
            this.sep = options.sep

        if( options.strip == true )
            this.stripBlanks = true

        // the header: can be a boolean or the list of columns
        if( options.header ) {
            if( options.header == true )
                firstLineAsHeader = true
            else if( options.header instanceof List )
                columnsHeader = options.header as List
            else
                throw new IllegalArgumentException("Not a valid header parameter value: ${options.header}")
        }

        // the quote character if used
        if( options.quote )
            quote = options.quote

        if( options.skip )
            skipLines = options.skip as int

        return this
    }

    /**
     * @return A map representing the valid options for the splitter. The map keys define the
     * accepted parameter names, the values the valid values for each of them.
     */
    protected Map<String,?> validOptions() {
        def result = super.validOptions()
        result.remove('file')       // <-- `file` mode not support by CsvSplitter
        result.remove('compress')   // <-- `compress` mode not supported
        result.sep = String
        result.strip = Boolean
        result.header = [ Boolean, List ]
        result.quote = String
        result.skip = Integer
        return result
    }

    /**
     * Implements the CSV parsing
     *
     * @param targetObject
     * @param index
     * @return
     */
    @Override
    protected process( Reader targetObject ) {
        final reader = wrapReader(targetObject)
        parseHeader(reader)
        super.process(reader)
    }

    protected void parseHeader(BufferedReader reader) {

        String line
        int z = 0

        while( z++ < skipLines && reader.readLine()) { /* nope */ }

        if( firstLineAsHeader ) {
            line = reader.readLine()
            if( !line ) throw new IllegalStateException("Missing 'header' in CSV file")
            String[] allCols = StringUtils.splitPreserveAllTokens(line, sep)
            columnsHeader = new ArrayList<>(allCols.size())
            for( int i=0; i<allCols.size(); i++ ) {
                def col = strip(allCols[i])
                if( !col ) throw new IllegalStateException("Empty header columns are not allowed in CSV file")
                columnsHeader[i] = col
            }
        }
    }

    /**
     * Process a CSV row at time
     *
     * @param reader The reader from which read the row
     * @return A {@link List} or a {@link Map} object representing the CSV row
     */
    @Override
    protected fetchRecord(BufferedReader reader) {
        String line
        while( true ) {
            line = reader.readLine()
            if( line )
                break
            if( line==null )
                return null
        }

        def tokens = StringUtils.splitPreserveAllTokens(line, sep) as List<String>
        // -- strip blanks and quote
        for( int i=0; i<tokens.size(); i++ ) {
            tokens[i] = strip(tokens[i])
        }

        // -- convert to a map if there's a columns header
        if( !columnsHeader )
            return tokens

        def map = [:]
        for( int i=0; i<columnsHeader.size(); i++ )
            map[ columnsHeader[i] ] = i<tokens.size() ? tokens[i] : null

        return map
    }

    /**
     * Remove the quote characters and extra blanks
     *
     * @param str The string to be stripped
     * @return The resulting string value
     */
    @PackageScope
    final String strip( String str ) {
        def value = stripBlanks ? StringUtils.strip(str) : str

        if( !quote )
            return value

        if( value.length()>1 && value.startsWith(quote) && value.endsWith(quote) )
            return value.substring(1, value.length()-1)

        else
            return value
    }

    protected CollectorStrategy createCollector() {

        if( counter.isEnabled() )
            return new ObjectListCollector()

        return null
    }

}
