package nextflow.splitter

import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.extension.FilesExtensions
/**
 * Splits a generic byte array in chunks having the specified length
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@InheritConstructors
class BytesSplitter extends AbstractBinarySplitter {

    @Override
    def apply( InputStream targetObject, int index ) {
        assert targetObject != null

        def result = null
        int c=0
        byte[] buffer = new byte[count]
        int item

        try {

            while( (item=targetObject.read()) != -1 ) {
                buffer[c] = (byte)item

                if ( ++c == count ) {
                    c = 0
                    result = invokeEachClosure(closure, buffer, index++ )
                    if( into != null )
                        append(into,result)

                    buffer = new byte[count]
                }
            }

        }
        finally {
            FilesExtensions.closeQuietly(targetObject)
        }


        /*
         * if there's something remaining in the buffer it's supposed
         * to be the last entry
         */

        if ( c ) {
            if( c != count ) {
                def copy = new byte[c]
                System.arraycopy(buffer,0,copy,0,c)
                buffer = copy
            }

            result = invokeEachClosure(closure, buffer, index++ )
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
            append(into,PoisonPill.instance)
            return into
        }
        if( into != null )
            return into

        return result
    }

}
