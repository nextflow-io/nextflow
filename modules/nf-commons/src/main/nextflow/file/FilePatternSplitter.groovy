/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.file
import java.nio.file.Path
import java.util.regex.Pattern

import groovy.transform.CompileStatic
/**
 * Parse a file path to isolate the parent, file-name and whenever it contains a glob|regex pattern
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class FilePatternSplitter {

    static final public FilePatternSplitter GLOB = glob()

    static enum Syntax { GLOB, REGEX }

    static final public Pattern GLOB_CURLY_BRACKETS = Pattern.compile(/(.*)(\{.*,.*\})(.*)/)

    static final public Pattern GLOB_SQUARE_BRACKETS = Pattern.compile(/(.*)(\[.+\])(.*)/)

    static final private char BACK_SLASH = '\\' as char

    static final private String GLOB_CHARS = '*?[]{}'

    static final private String REGEX_CHARS = '.^$+{}[]|()'

    private boolean pattern

    private final Syntax syntax

    private String scheme

    private String parent

    private String fileName

    boolean isPattern() { pattern }

    String getFileName() { fileName }

    String getParent() { parent }

    String getScheme() { scheme }

    static FilePatternSplitter glob() { new FilePatternSplitter(Syntax.GLOB) }

    static FilePatternSplitter regex() { new FilePatternSplitter(Syntax.REGEX) }

    FilePatternSplitter( Syntax syntax ) {
        this.syntax = syntax
    }

    private String metaChars() {
        syntax == Syntax.GLOB ? GLOB_CHARS : REGEX_CHARS
    }

    private boolean containsMetaChars(String str) {
        final meta = metaChars()

        for( int i=0; i<str.length(); i++ ) {
            if( meta.contains(str[i]) )
                return true
        }

        return false
    }

    /**
     * Parse a file path detecting the schema, parent folder, file name and pattern components
     *
     * @param filePath The file path string to parse
     * @return
     */
    FilePatternSplitter parse(String filePath) {
        // -- detected the file scheme if any
        int p = filePath.indexOf('://')
        if( p != -1 ) {
            scheme = filePath.substring(0, p)
            filePath = filePath.substring(p+3)
        }
        else {
            scheme = null
        }

        //
        // split the path in two components
        // - folder: the part not containing meta characters
        // - pattern: the part containing a meta character

        boolean found
        String norm = replaceMetaChars(filePath)
        p = firstMetaIndex(norm)

        if( p == -1 ) {
            // find the last SLASH
            p = filePath.lastIndexOf('/')
            found = false
        }
        else {
            found = true
            // walk back to the first SLASH char
            int i = p
            p = -1
            while( --i >= 0 ) {
                if( filePath[i] == '/' )  {
                    p = i
                    break
                }
            }
        }

        if( p == -1 ) {
            parent = './'
            fileName = filePath
            pattern = found && pairedBrackets(norm)
        }
        else {
            parent = strip(filePath.substring(0,p+1))
            fileName = filePath.substring(p+1)
            pattern = found && pairedBrackets(norm)
        }

        return this
    }


    private boolean pairedBrackets(String str) {
        if( syntax == Syntax.REGEX )
            return true

        if( str.contains('{') || str.contains('}') )
            return GLOB_CURLY_BRACKETS.matcher(str).matches()

        if( str.contains('[') || str.contains(']') )
            return GLOB_SQUARE_BRACKETS.matcher(str).matches()

        return true
    }

    protected String replaceMetaChars( String str, char marker = 0x0 ) {
        // create a version replacing escape meta chars with an 0x0
        final meta = metaChars()
        final result = new StringBuilder()
        int i=0;
        while( i<str.length() ) {
            def ch = str.charAt(i++)
            if( ch == BACK_SLASH && i<str.length() && meta.contains(str[i])) {
                result.append(ch).append(marker)
                i++
            }
            else {
                result.append(ch)
            }
        }

        return result
    }

    protected int firstMetaIndex(String str) {

        // find the index of the first meta chars
        final meta = metaChars()
        def min = Integer.MAX_VALUE
        for( int i=0; i<meta.length(); i++ ) {
            def p = str.indexOf(meta[i])
            if( p!=-1 && p<min )
                min = p
        }

        return min != Integer.MAX_VALUE ? min : -1
    }


    /**
     * Strips backslash characters from a path
     */
    String strip( String str ) {
        int p = str.indexOf('\\')
        if( p == -1 )
            return str

        final meta = metaChars()
        final result = new StringBuilder()
        int i=0;
        while( i<str.length() ) {
            def ch = str.charAt(i++)
            if( ch != BACK_SLASH || i==str.length() || !meta.contains(str[i])) {
                result.append(ch)
            }
        }

        return result.toString()
    }

    /**
     * Escape pattern meta-characters with a backslash
     */
    String escape( Path path ) {
        escape(path.toString())
    }

    /**
     * Escape pattern meta-characters with a backslash
     */
    String escape( String str ) {
        final meta = metaChars()
        def result = str
        for( int i=0; i<meta.length(); i++ ) {
            result = result.replace( meta[i], '\\' + meta[i] )
        }
        return result
    }
}
