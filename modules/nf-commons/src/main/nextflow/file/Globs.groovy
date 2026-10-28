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

package nextflow.file

import java.util.regex.PatternSyntaxException

import groovy.transform.CompileStatic

/**
 * Convert a glob pattern into a regular expression, following the
 * glob syntax of {@link java.nio.file.FileSystem#getPathMatcher(java.lang.String)}
 * on Unix file systems
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Globs {

    /**
     * Translate a glob pattern into a regular expression that
     * must match the entire path string
     *
     * @param glob The glob pattern
     * @return The equivalent regular expression
     * @throws PatternSyntaxException if the glob pattern is invalid
     */
    static String toUnixRegexPattern(String glob) {
        return new GlobTranslator(glob).translate()
    }

    private static class GlobTranslator {

        private static final String REGEX_SPECIAL = '\\^$.|?*+()[]{}'

        private static final String CLASS_SPECIAL = '\\[]&-^'

        private final String glob

        private final StringBuilder regex = new StringBuilder()

        private int pos

        GlobTranslator(String glob) {
            this.glob = glob
        }

        String translate() {
            regex.append('^')
            while( !atEnd() )
                translateElement()
            regex.append('$')
            return regex.toString()
        }

        private boolean atEnd() {
            return pos >= glob.length()
        }

        private char peek(int offset=0) {
            return glob.charAt(pos + offset)
        }

        private void translateElement() {
            final char ch = peek()
            if( ch == '*' as char )
                translateWildcard()
            else if( ch == '?' as char ) {
                pos++
                regex.append('[^/]')
            }
            else if( ch == '[' as char )
                translateClass()
            else if( ch == '{' as char )
                translateGroup()
            else if( ch == '\\' as char )
                translateEscape()
            else {
                pos++
                appendLiteral(ch)
            }
        }

        private void translateWildcard() {
            if( pos + 1 < glob.length() && peek(1) == '*' as char ) {
                // `**` crosses directory boundaries but, as with the JDK, not line terminators
                pos += 2
                regex.append('.*')
            }
            else {
                pos++
                regex.append('[^/]*')
            }
        }

        private void translateEscape() {
            final int start = pos++
            if( atEnd() )
                throw syntaxError('Missing character after escape', start)
            appendLiteral(glob.charAt(pos++))
        }

        private void translateGroup() {
            final int start = pos++
            regex.append('(?:')
            while( true ) {
                if( atEnd() )
                    throw syntaxError('Unclosed group', start)
                final char ch = peek()
                if( ch == '}' as char ) {
                    pos++
                    break
                }
                if( ch == ',' as char ) {
                    pos++
                    regex.append('|')
                }
                else if( ch == '{' as char )
                    throw syntaxError('Nested groups are not supported', pos)
                else
                    translateElement()
            }
            regex.append(')')
        }

        private void translateClass() {
            final int start = pos++
            final boolean negated = !atEnd() && peek() == '!' as char
            if( negated )
                pos++
            final int first = pos
            final StringBuilder members = new StringBuilder()
            boolean spansSlash = false
            // last single character, which may be the start of a range, or -1
            int rangeStart = -1

            while( true ) {
                if( atEnd() )
                    throw syntaxError('Unclosed character class', start)
                final char ch = peek()
                if( ch == ']' as char ) {
                    if( pos == first )
                        throw syntaxError('Empty character class', start)
                    pos++
                    break
                }
                if( ch == '/' as char )
                    throw syntaxError('Character class cannot match a path separator', pos)
                if( ch == '-' as char && pos != first ) {
                    if( rangeStart < 0 )
                        throw syntaxError('Invalid character range', pos)
                    if( pos + 1 >= glob.length() )
                        throw syntaxError('Unclosed character class', start)
                    final char end = peek(1)
                    if( end == ']' as char ) {
                        // a trailing `-` is literal
                        members.append('\\-')
                        pos++
                        continue
                    }
                    // the JDK also rejects `\` and `[` as the end of a range
                    if( end < rangeStart || end == '\\' as char || end == '[' as char )
                        throw syntaxError('Invalid character range', pos)
                    if( rangeStart < ('/' as char) && end >= ('/' as char) )
                        spansSlash = true
                    members.append('-')
                    appendClassMember(members, end)
                    rangeStart = -1
                    pos += 2
                    continue
                }
                appendClassMember(members, ch)
                // like a leading `-`, a leading `^` cannot start a range (as with the JDK)
                final boolean leading = pos == first && (ch == '-' as char || (ch == '^' as char && !negated))
                rangeStart = leading ? -1 : (int) ch
                pos++
            }

            if( negated )
                regex.append('[^/').append(members).append(']')
            else if( spansSlash )
                regex.append('(?!/)[').append(members).append(']')
            else
                regex.append('[').append(members).append(']')
        }

        private void appendLiteral(char ch) {
            if( REGEX_SPECIAL.indexOf((int) ch) >= 0 )
                regex.append('\\')
            regex.append(ch)
        }

        private static void appendClassMember(StringBuilder members, char ch) {
            if( CLASS_SPECIAL.indexOf((int) ch) >= 0 )
                members.append('\\')
            members.append(ch)
        }

        private PatternSyntaxException syntaxError(String message, int index) {
            return new PatternSyntaxException(message, glob, index)
        }
    }

}
