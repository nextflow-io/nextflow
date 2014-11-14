package nextflow.splitter
import java.nio.charset.Charset
import java.nio.charset.CharsetDecoder
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import nextflow.util.CharsetHelper
/**
 * Implements an abstract text splitting for the main types
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@CompileStatic
@InheritConstructors
abstract class AbstractTextSplitter extends AbstractSplitter<Reader> {

    protected Charset charset = Charset.defaultCharset()

    Charset getCharset() { charset }

    AbstractTextSplitter options(Map options) {
        super.options(options)
        charset = CharsetHelper.getCharset(options.charset)
        return this
    }

    /**
     * @return A map representing the valid options for the splitter. The map keys define the
     * accepted parameter names, the values the valid values for each of them.
     */
    @Override
    protected Map<String,?> validOptions() {
        def result = super.validOptions()
        result.charset = [ Charset, Map, String ]
        return result
    }


    protected Reader normalizeType( obj ) {

        if( obj instanceof Reader )
            return (Reader) obj

        if( obj instanceof CharSequence )
            return new StringReader(obj.toString())

        if( obj instanceof Path )
            return newReader(obj, charset)

        if( obj instanceof InputStream )
            return new InputStreamReader(obj,charset)

        if( obj instanceof File )
            return newReader(obj.toPath(), charset)

        if( obj instanceof char[] )
            return new StringReader(new String((char[])obj))

        throw new IllegalArgumentException("Object of class '${obj.class.name}' does not support 'splitter' methods")

    }

    protected Reader newReader( Path path, Charset charset ) {
        def source = newInputStream(path)
        CharsetDecoder decoder = charset.newDecoder()
        Reader reader = new InputStreamReader(source, decoder)
        new BufferedReader(reader)
    }


}
