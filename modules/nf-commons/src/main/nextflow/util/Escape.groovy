/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.util

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.extension.FilesEx

/**
 * Escape helper class
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Escape {

    final private static List<String> SPECIAL_CHARS = ["'", '"', ' ', '(', ')', '\\', '!', '&', '|', '<', '>', '`', ':', ';']

    final private static List<String> VAR_CHARS = ['$']

    // a task data file name is literal, so `$` is escaped too
    final private static List<String> DATA_PATH_CHARS = SPECIAL_CHARS + ['$']

    final private static List<String> WILDCARDS = ["*", "?", "{", "}", "[", "]", "'", '"', ' ', '(', ')', '\\', '!', '&', '|', '<', '>', '`', ':']

    private static String replace(List<String> special, String str, boolean doNotEscapeComplement=false) {
        def copy = new StringBuilder(str.size() +10)
        for( int i=0; i<str.size(); i++) {
            def ch = str[i]
            def p = special.indexOf(ch)
            if( p != -1 ) {
                // when ! is the first character after a `[` it should not be escaped
                // see http://man7.org/linux/man-pages/man7/glob.7.html
                final isComplement = doNotEscapeComplement && ch=='!' && ( i>0 && str[i-1]=='[' && (i==1 || str[i-2]!='\\') && str.substring(i).contains(']'))
                if( !isComplement )
                    copy.append('\\')
            }
            copy.append(str[i])
        }
        return copy.toString()
    }

    static String wildcards(String str) {
        replace(WILDCARDS, str)
    }

    static String path(String val) {
        replace(SPECIAL_CHARS, val, true)
    }

    static String path(Path val) {
        path(val.toString())
    }

    /**
     * Escape a path that comes from the task data, e.g. an input file name, so that it is used literally
     */
    static String dataPath(String val) {
        checkNewline(val)
        replace(DATA_PATH_CHARS, val, true)
    }

    static String dataPath(Path val) {
        dataPath(val.toString())
    }

    /**
     * Escape a task data file name used in a task script, where `$` is kept as is so the name also
     * works in quoted strings and non-Bash scripts
     */
    static String scriptPath(String val) {
        checkNewline(val)
        // `${...}` and `$[...]` are expanded by Bash even when the name is escaped by `path`, so they are not supported
        if( val.contains('${') || val.contains('$[') )
            throw new IllegalArgumentException("File name cannot contain '\${' or '\$[': $val")
        path(val)
    }

    // a newline cannot be kept in a single shell word in every context, so such names are not supported
    private static void checkNewline(String val) {
        if( val.contains('\n') )
            throw new IllegalArgumentException("File name cannot contain a newline character: ${val.replace('\n','\\n')}")
    }

    static String path(File val) {
        path(val.toString())
    }

    static String path(GString val) {
        path(val.toString())
    }

    static String cli(String[] args) {
        args.collect { cli(it) }.join(' ')
    }

    static String cli(String arg) {
        if( arg.contains("'") )
            return wildcards(arg)
        def x = wildcards(arg)
        x == arg ? arg : "'" + arg + "'"
    }

    static String uriPath(Path file) {
        final scheme = FilesEx.getScheme(file)
        if( scheme == 'file' )
            return path(file)

        final uri = FilesEx.toUriString(file)
        final prefix = "$scheme://"
        assert uri.startsWith(prefix)
        def loc = uri.substring(prefix.length())
        while( loc.size()>1 && loc.endsWith('/') )
            loc = loc.substring(0,loc.length()-1)
        return prefix + path(loc)
    }

    static String blanks(String str) {
        str
                .replaceAll('\n',/\\n/)
                .replaceAll('\t',/\\t/)
                .replaceAll('\r',/\\r/)
                .replaceAll('\f',/\\f/)

    }

    static String variable(String val) {
        replace(VAR_CHARS, val, false)
    }

    /**
     * Escape the characters that are special inside a bash double-quoted string, e.g. the argument of `eval "..."`
     */
    static String doubleQuoted(String val) {
        val.replaceAll('[\\\\"$`]', '\\\\$0')
    }
}
