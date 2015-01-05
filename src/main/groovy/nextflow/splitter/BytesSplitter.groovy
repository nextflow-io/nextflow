package nextflow.splitter
import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
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
    def process( InputStream targetObject) {
        assert targetObject != null

        def result = null
        int c=0
        long bytesCount=0
        byte[] buffer = new byte[count]
        int item

        try {

            while( (item=targetObject.read()) != -1 ) {
                buffer[c] = (byte)item

                if ( ++c == count ) {
                    c = 0
                    result = invokeEachClosure(closure, buffer)
                    buffer = new byte[count]
                }

                // -- check the limit of allowed rows has been reached
                if( limit && ++bytesCount == limit )
                    break
            }

        }
        finally {
            targetObject.closeQuietly()
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

            result = invokeEachClosure(closure, buffer)
        }

        return result
    }

    @Override
    protected CollectorStrategy createCollector() {
        return null
    }
}
