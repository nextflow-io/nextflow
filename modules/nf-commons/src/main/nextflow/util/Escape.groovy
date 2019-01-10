/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

/**
 * Escape helper class
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Escape {

    private static List<String> SPECIAL_CHARS = ["'", '"', ' ', '(', ')', '\\', '!', '&', '|', '<', '>', '`', ':']

    private static List<String> WILDCARDS = ["*", "?", "{", "}", "[", "]", "'", '"', ' ', '(', ')', '\\', '!', '&', '|', '<', '>', '`', ':']

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

}
