package nextflow.splitter
import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.PoisonPill
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

    @Override
    def apply( Reader targetObject, int index ) {

        def result = null
        def buffer = new StringBuilder()
        int c = 0
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
