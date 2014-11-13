package nextflow.splitter
import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.PoisonPill
/**
 * Split FASTQ formatted text content or files
 *
 * @link http://en.wikipedia.org/wiki/FASTQ_format
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@InheritConstructors
class FastqSplitter extends AbstractTextSplitter {

    static Map recordToMap( String l1, String l2, String l3, String l4, map ) {
        def result = [:]

        final isMap = map instanceof Map
        if( !isMap || (map as Map).containsKey('readHeader'))
            result.readHeader = l1.substring(1)

        if( !isMap || (map as Map).containsKey('readString'))
            result.readString = l2

        if( !isMap || (map as Map).containsKey('qualityHeader'))
            result.qualityHeader = l3.substring(1)

        if( !isMap || (map as Map).containsKey('qualityString'))
            result.qualityString = l4

        return result
    }

    static void recordToText( String l1, String l2, String l3, String l4, StringBuilder buffer ) {
        // read header
        buffer << l1 << '\n'
        // read string
        buffer << l2 << '\n'
        // quality header
        buffer << l3 << '\n'
        // quality string
        buffer << l4 << '\n'
    }

    @Override
    def apply( Reader targetObject, int index )  {

        BufferedReader reader0 = (BufferedReader)(targetObject instanceof BufferedReader ? targetObject : new BufferedReader(targetObject))

        final StringBuilder buffer = new StringBuilder()
        int blockCount=0
        def result = null

        def error = "Invalid FASTQ format"
        if( sourceFile )
            error += " for file: " + sourceFile

        try {
            while( true ) {
                def l1 = reader0.readLine()
                def l2 = reader0.readLine()
                def l3 = reader0.readLine()
                def l4 = reader0.readLine()

                if( !l1 || !l2 || !l3 || !l4 )
                    break

                if( !l1.startsWith('@') || !l3.startsWith('+') )
                    throw new IllegalStateException(error)

                if( !recordMode )
                    recordToText(l1,l2,l3,l4,buffer)

                if ( ++blockCount == count ) {
                    // invoke the closure, passing the read block as parameter
                    def splitArg = recordMode ? recordToMap(l1,l2,l3,l4, recordCols) : buffer.toString()

                    result = invokeEachClosure( closure, splitArg, index++ )
                    if( into != null ) {
                        append(into,result)
                    }

                    buffer.setLength(0)
                    blockCount=0
                }

            }
        }
        finally {
            reader0.closeQuietly()
        }

        /*
         * if there's something remaining in the buffer it's supposed
         * to be the last entry
         */
        if ( buffer.size() ) {
            // invoke the closure, passing the read block as parameter
            result = invokeEachClosure( closure, buffer.toString(), index )
            if( into != null )
                append(into,result)
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
}
