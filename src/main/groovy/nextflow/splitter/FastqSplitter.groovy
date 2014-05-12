package nextflow.splitter

import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.PoisonPill
import net.sf.picard.fastq.FastqReader
import net.sf.picard.fastq.FastqRecord
import nextflow.extension.FilesExtensions
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@InheritConstructors
class FastqSplitter extends AbstractTextSplitter {

    static Map recordToMap( FastqRecord record, map ) {
        def result = [:]

        final isMap = map instanceof Map
        if( !isMap || (map as Map).containsKey('readHeader'))
            result.readHeader = record.readHeader

        if( !isMap || (map as Map).containsKey('readString'))
            result.readString = record.readString

        if( !isMap || (map as Map).containsKey('baseQualityHeader'))
            result.qualityHeader = record.baseQualityHeader

        if( !isMap || (map as Map).containsKey('baseQualityString'))
            result.qualityString = record.baseQualityString

        return result
    }

    static void recordToText( FastqRecord record, StringBuilder buffer ) {
        // read header
        buffer << '@' << record.readHeader << '\n'
        // read string
        buffer << record.readString << '\n'
        // quality header
        buffer << '+'
        if( record.baseQualityHeader ) buffer << record.baseQualityHeader
        buffer << '\n'
        // quality string
        buffer << record.baseQualityString
        buffer << '\n'
    }

    @Override
    def apply( Reader targetObject, int index )  {

        BufferedReader reader0 = (BufferedReader)(targetObject instanceof BufferedReader ? targetObject : new BufferedReader(targetObject))

        final fastq = new FastqReader(reader0)
        final StringBuilder buffer = new StringBuilder()
        int blockCount=0
        def result = null

        try {
            while( fastq.hasNext() ) {
                def record = fastq.next()

                if( !recordMode )
                    recordToText(record,buffer)

                if ( ++blockCount == count ) {
                    // invoke the closure, passing the read block as parameter
                    def splitArg = recordMode ? recordToMap(record, recordCols) : buffer.toString()

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
            FilesExtensions.closeQuietly(fastq)
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
