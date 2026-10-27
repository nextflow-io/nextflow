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
import nextflow.file.FilePatternSplitter

/**
 * Escape helper class
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Escape {

    final private static List<String> SPECIAL_CHARS = ["'", '"', ' ', '(', ')', '\\', '!', '&', '|', '<', '>', '`', ':', ';', '\n']

    final private static List<String> VAR_CHARS = ['$']

    final private static List<String> WILDCARDS = ["*", "?", "{", "}", "[", "]", "'", '"', ' ', '(', ')', '\\', '!', '&', '|', '<', '>', '`', ':', '\n']

    private static String replace(List<String> special, String str, boolean doNotEscapeComplement=false, boolean keepEscapes=false) {
        def copy = new StringBuilder(str.size() +10)
        for( int i=0; i<str.size(); i++) {
            def ch = str[i]
            // in a glob pattern a backslash makes the next character literal, so keep the pair
            if( keepEscapes && ch == '\\' && i+1 < str.size() ) {
                final next = str[++i]
                copy.append( next == '\n' ? '$\'\\n\'' : '\\' + next )
                continue
            }
            def p = special.indexOf(ch)
            // a backslash-newline is a line continuation, so use ANSI-C quoting; it also survives `eval "..."`
            if( p != -1 && ch == '\n' ) {
                copy.append('$\'\\n\'')
                continue
            }
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

    /**
     * Escape an output file pattern the way the output file collector reads it:
     * a glob pattern keeps its backslash escapes, e.g. `file\[1\].txt`, any other value is a literal file name
     */
    static String pattern(String val) {
        final splitter = FilePatternSplitter.glob().parse(val)
        return splitter.isPattern()
            ? replace(SPECIAL_CHARS, val, true, true)
            : path(splitter.strip(val))
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
