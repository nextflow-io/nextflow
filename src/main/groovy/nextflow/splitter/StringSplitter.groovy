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

    protected CollectorStrategy createCollector() {
        count > 1 ? new CharSequenceCollector() : null
    }

    @Override
    protected fetchRecord(BufferedReader targetObject) {

        int ch
        while( true ) {
            ch = targetObject.read()
            if( ch == -1 )
                return null

            if( ignoreNewLine && ( ch == '\n' as char || ch == '\r' as char ))
                continue

            break
        }

        return ch as char
    }
}
