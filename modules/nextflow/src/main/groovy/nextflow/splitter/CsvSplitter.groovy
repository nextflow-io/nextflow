/*
 * Copyright 2013-2023, Seqera Labs
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
import groovy.util.logging.Slf4j
import nextflow.util.CsvParser
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
     * Implement the CSV splitting logic
     */
    private CsvParser parser = new CsvParser(empty: '')

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
     * Define the types to cast each csv column.
     */
    static private final Set<String> VALID_COLUMN_TYPES = ['string', 'boolean', 'character', 'short', 'integer', 'long', 'float', 'double'] as Set

    protected List<String> columnTypes
    protected Map<String,String> columnTypesMap
    protected Closure castMethod

    /**
     * Set the splitter options by specifying a map of named parameters.
     * Valid parameters are:
     * <li>{@code sep}
     * <li>{@code strip}
     * <li>{@code header}
     * <li>{@code quote}
     * <li>{@code skip}
     * <li>{@code types}
     * <li>{@code strict}
     *
     * @param options
     * @return The splitter instance itself
     */
    CsvSplitter options(Map options) {

        super.options(options)

        // the separator character
        if( options.sep )
            parser.setSeparator(options.sep as String)

        if( options.strip == true )
            parser.setStrip(true)

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
            parser.setQuote(options.quote as String)

        if( options.skip )
            skipLines = options.skip as int

        // cast variables to their type
        castMethod = options.strict
            ? this.&castValueTypeStrict
            : this.&castValueType

        if( options.types instanceof List ) {
            if( !options.types.every { VALID_COLUMN_TYPES.contains(it) } )
                throw new IllegalArgumentException("Provided types are not allowed: ${options.types}. Valid column types are: ${VALID_COLUMN_TYPES}")

            columnTypes = options.types as List
        }
        else if( options.types instanceof Map ) {
            if( !options.header )
                throw new IllegalArgumentException("Column types map can only be used with the `header` option.")

            if( !options.types.every { VALID_COLUMN_TYPES.contains(it.value) } )
                throw new IllegalArgumentException("Provided types are not allowed: ${options.types}. Valid column types are: ${VALID_COLUMN_TYPES}")

            columnTypesMap = options.types as Map
        }
        else if( options.types != null )
            throw new IllegalArgumentException("Not a valid types parameter value: ${options.types}")

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
        result.types = [ List, Map ]
        result.strict = Boolean
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
        // parse header row if specified
        final reader = wrapReader(targetObject)
        parseHeader(reader)

        // convert columns map to columns list
        if( columnsHeader && columnTypesMap )
            columnTypes = columnsHeader.collect { columnTypesMap[it] ?: 'string' }

        super.process(reader)
    }

    protected void parseHeader(BufferedReader reader) {

        String line
        int z = 0

        while( z++ < skipLines && reader.readLine()) { /* nope */ }

        if( firstLineAsHeader ) {
            line = reader.readLine()
            if( !line ) throw new IllegalStateException("Missing 'header' in CSV file")
            List allCols = parser.parse(line)
            columnsHeader = new ArrayList<>(allCols.size())
            for( int i=0; i<allCols.size(); i++ ) {
                def col = allCols[i]
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

        final tokens = parser.parse(line)

        // -- convert to a map if there's a columns header
        if( !columnsHeader )
            return tokens

        if( columnTypes && columnTypes.size() != tokens.size() )
            throw new IllegalArgumentException("The number of column types should match the number of csv columns. Provided types: ${columnTypes}")

        def map = [:]
        for( int i = 0; i < tokens.size(); i++ )
            map[columnsHeader[i]] = columnTypes
                ? castMethod(tokens[i], columnTypes[i])
                : tokens[i]

        for( int i = tokens.size(); i < columnsHeader.size(); i++ )
            map[columnsHeader[i]] = null

        return map
    }

    /**
     * Cast a value to the provided type
     *
     * @param str
     * @param type
     * @return The value casted to its primitive type
     */
    static protected castValueType(String str, String type) {

        try {
            if( str == null || str == '' ) return null
            else if( type.toLowerCase() == 'boolean' ) return str.toBoolean()
            else if( type.toLowerCase() == 'character' ) return str.toCharacter()
            else if( type.toLowerCase() == 'short' ) return str.toShort()
            else if( type.toLowerCase() == 'integer' ) return str.toInteger()
            else if( type.toLowerCase() == 'long' ) return str.toLong()
            else if( type.toLowerCase() == 'float' ) return str.toFloat()
            else if( type.toLowerCase() == 'double' ) return str.toDouble()
            else if( type.toLowerCase() == 'string' ) return str
            else {
                log.warn "Unable to cast value $str to type $type"
                return str
            }
        } catch( Exception e ) {
            log.warn "Unable to cast value $str to type $type: $e"
            return str
        }
    }

    /**
     * Cast a value to the provided type in a Strict mode
     *
     * @param str
     * @param type
     * @return The value casted to its primitive type
     */
    static private final Set<String> VALID_BOOLEAN_VALUES = ['true', 'false'] as Set

    static protected castValueTypeStrict(String str, String type) {

        try {
            if( str == null || str == '' ) return null
            else if( type.toLowerCase() == 'boolean' && str.toLowerCase() in VALID_BOOLEAN_VALUES ) return str.toBoolean()
            else if( type.toLowerCase() == 'character' ) return str.toCharacter()
            else if( type.toLowerCase() == 'short' && str.isNumber() ) return str.toShort()
            else if( type.toLowerCase() == 'integer' && str.isInteger() ) return str.toInteger()
            else if( type.toLowerCase() == 'long' && str.isLong() ) return str.toLong()
            else if( type.toLowerCase() == 'float' && str.isFloat() ) return str.toFloat()
            else if( type.toLowerCase() == 'double' && str.isDouble() ) return str.toDouble()
            else if( type.toLowerCase() == 'string' ) return str
            else {
                log.warn "Value $str is not of type $type: returning a string"
                return str
            }
        } catch( Exception e ) {
            log.warn "Unable to cast value $str to type $type: $e"
            return str
        }
    }

    protected CollectorStrategy createCollector() {

        if( counter.isEnabled() )
            return new ObjectListCollector()

        return null
    }

}
