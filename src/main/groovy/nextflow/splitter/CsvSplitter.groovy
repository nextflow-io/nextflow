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
    def process( Reader targetObject, int index ) {
        BufferedReader reader0 = (BufferedReader)(targetObject instanceof BufferedReader ? targetObject : new BufferedReader(targetObject))

        def result = null
        String line
        int c=0
        long itemsCount=0

        List buffer = count > 1 ? [] : null

        int z = 0
        while( z++ < skipLines && reader0.readLine()) { /* nope */ }


        if( firstLineAsHeader ) {
            line = reader0.readLine()
            if( !line ) throw new IllegalStateException("Missing 'header' in CSV file")
            String[] allCols = StringUtils.splitPreserveAllTokens(line, sep)
            columnsHeader = new ArrayList<>(allCols.size())
            for( int i=0; i<allCols.size(); i++ ) {
                def col = strip(allCols[i])
                if( !col ) throw new IllegalStateException("Empty header columns are not allowed in CSV file")
                columnsHeader[i] = col
            }
        }

        try {
            while( (line = reader0.readLine()) != null ) {

                def row
                def tokens = StringUtils.splitPreserveAllTokens(line, sep)
                // -- strip blanks and quote
                for( int i=0; i<tokens.length; i++ ) {
                    tokens[i] = strip(tokens[i])
                }

                // -- convert to a map if there's a columns header
                if( !columnsHeader ) {
                    row = tokens
                }
                else {
                    def map = [:]
                    for( int i=0; i<columnsHeader.size(); i++ )
                        map[ columnsHeader[i] ] = i<tokens.size() ? tokens[i] : null
                    row = map
                }


                // -- append to the list buffer
                if( buffer != null ) {
                    buffer.add(row)
                }

                // -- invoke the closure
                if ( ++c == count ) {
                    c = 0
                    result = invokeEachClosure(closure, buffer ?: row, index++ )
                    if( into != null )
                        append(into,result)

                    if( buffer!=null )
                        buffer = []
                }

                // -- check the limit of allowed rows has been reached
                if( limit && ++itemsCount == limit )
                    break
            }

        }
        finally {
            reader0.closeQuietly()
        }

        /*
         * if there's something remaining in the buffer it's supposed
         * to be the last entry
         */
        if ( buffer ) {
            result = invokeEachClosure(closure, buffer, index )
            if( into != null )
                append(into, result)
        }

        return result
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

}
