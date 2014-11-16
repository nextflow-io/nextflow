package nextflow.splitter
import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
/**
 * Split a text file by one or more lines at times
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@InheritConstructors
class TextSplitter extends AbstractTextSplitter {

    @Override
    def process( Reader targetObject, int index )  {
        BufferedReader reader0 = (BufferedReader)(targetObject instanceof BufferedReader ? targetObject : new BufferedReader(targetObject))

        def result = null
        String line
        StringBuilder buffer = new StringBuilder()
        int c=0
        long itemsCount=0

        try {

            while( (line = reader0.readLine()) != null ) {
                buffer << line << '\n'

                if ( ++c == count ) {
                    c = 0
                    result = invokeEachClosure(closure, buffer.toString(), index++ )
                    if( into != null )
                        append(into,result)

                    buffer.setLength(0)
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
            result = invokeEachClosure(closure, buffer.toString(), index )
            if( into != null )
                append(into, result)
        }

        return result
    }

}
