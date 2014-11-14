package nextflow.splitter
import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
/**
 * Simple slitter chunking a string in sub-strings having the specified length
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@InheritConstructors
class StringSplitter extends AbstractTextSplitter {

    protected boolean ignoreNewLine

    StringSplitter options(Map options) {
        super.options(options)
        ignoreNewLine = options.ignoreNewLine == true ?: false
        return this
    }

    /**
     * @return A map representing the valid options for the splitter. The map keys define the
     * accepted parameter names, the values the valid values for each of them.
     */
    @Override
    protected Map<String,?> validOptions() {
        def result = super.validOptions()
        result.ignoreNewLine = Boolean
        return result
    }

    @Override
    def process( Reader targetObject, int index ) {

        def result = null
        def buffer = new StringBuilder()
        int c = 0
        long itemsCount=0
        def ch

        try {

            while( (ch=targetObject.read()) != -1 ) {
                if( ignoreNewLine && ( ch == '\n' as char || ch == '\r' as char ))
                    continue
                buffer.append( (char)ch )

                if ( ++c == count ) {
                    c = 0
                    result = invokeEachClosure(closure, buffer.toString(), index++ )
                    if( into != null )
                        append(into, result)

                    buffer.setLength(0)
                }

                // -- check the limit of allowed rows has been reached
                if( limit && ++itemsCount == limit )
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
        if ( buffer.size() ) {
            result = invokeEachClosure(closure, buffer.toString(), index++ )
            if( into != null )
                append(into, result)
        }

        return result
    }
}
