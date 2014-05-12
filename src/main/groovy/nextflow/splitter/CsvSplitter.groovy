package nextflow.splitter

import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.extension.FilesExtensions
import org.apache.commons.lang.StringUtils

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
@InheritConstructors
class CsvSplitter extends AbstractTextSplitter {

    protected String sep = ','

    protected String quote

    protected boolean stripBlanks

    protected boolean firstLineAsHeader

    protected List<String> columnsHeader

    protected int skipLines = 0

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

    @Override
    def apply( Reader targetObject, int index ) {
        BufferedReader reader0 = (BufferedReader)(targetObject instanceof BufferedReader ? targetObject : new BufferedReader(targetObject))

        def result = null
        String line
        int c=0

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
                if( buffer!=null ) {
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
            }

        }
        finally {
            FilesExtensions.closeQuietly(reader0)
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

        /*
         * now close and return the result
         * - when the target it's a channel, send stop message
         * - when it's a list return it
         * - otherwise return the last value
         */
        if( into instanceof DataflowWriteChannel && autoClose ) {
            append(into, PoisonPill.instance)
            return into
        }
        if( into != null )
            return into

        return result
    }

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
