package nextflow.splitter
import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
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

    private boolean processQualityField

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
    def process( Reader targetObject, int index )  {

        BufferedReader reader0 = (BufferedReader)(targetObject instanceof BufferedReader ? targetObject : new BufferedReader(targetObject))

        final StringBuilder buffer = new StringBuilder()
        int blockCount=0
        long itemsCount=0
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
                    def closureArg
                    if( recordMode )
                        closureArg = recordToMap(l1,l2,l3,l4, recordFields)
                    else if( processQualityField )
                        closureArg = l4
                    else
                        closureArg = buffer.toString()

                    result = invokeEachClosure( closure, closureArg, index++ )
                    if( into != null ) {
                        append(into,result)
                    }

                    buffer.setLength(0)
                    blockCount=0
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
        if ( buffer.size() ) {
            // invoke the closure, passing the read block as parameter
            result = invokeEachClosure( closure, buffer.toString(), index )
            if( into != null )
                append(into,result)
        }

        return result
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
