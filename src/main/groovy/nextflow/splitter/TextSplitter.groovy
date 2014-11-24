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

    /**
     * A record is a text line
     *
     * @param reader The buffered reader
     * @return A line string or {@code null} when the end of the file is reached
     */
    @Override
    protected fetchRecord(BufferedReader reader) {
        def line = reader.readLine()
        if( line != null ) line+='\n'
        return line
    }
}
