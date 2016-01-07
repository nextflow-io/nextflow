package nextflow.util

import java.nio.file.Path

import groovy.transform.CompileStatic

/**
 * Escape helper class
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Escape {

    private static String QUOTE = "'"

    private static String DOUBLE_QUOTE = '"'

    private static String BLANK = " "

    static String path(String val) {
        val.replace(QUOTE,'\\'+QUOTE).replace(BLANK,'\\'+BLANK).replace(DOUBLE_QUOTE, '\\'+DOUBLE_QUOTE)
    }

    static String path(Path val) {
        path(val.toString())
    }

    static String path(File val) {
        path(val.toString())
    }

    static String path(GString val) {
        path(val.toString())
    }

}
