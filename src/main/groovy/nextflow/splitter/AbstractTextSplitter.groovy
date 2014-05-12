package nextflow.splitter
import java.nio.charset.Charset
import java.nio.file.Files
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

    protected Charset charset

    Charset getCharset() { charset }

    AbstractTextSplitter options(Map options) {
        super.options(options)
        charset = CharsetHelper.getCharset(options.charset)
        return this
    }

    protected Reader normalizeType( obj ) {

        if( obj instanceof Reader )
            return (Reader) obj

        if( obj instanceof CharSequence )
            return new StringReader(obj.toString())

        if( obj instanceof Path )
            return Files.newBufferedReader(obj, charset)

        if( obj instanceof InputStream )
            return new InputStreamReader(obj,charset)

        if( obj instanceof File )
            return new FileReader(obj)

        if( obj instanceof char[] )
            return new StringReader(new String(obj))

        throw new IllegalArgumentException("Object of class '${obj.class.name}' does not support 'chopString' method")

    }


}
