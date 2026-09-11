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

package nextflow.script.parser;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.antlr.v4.runtime.ParserRuleContext;
import org.codehaus.groovy.syntax.SyntaxException;

/**
 * Utilities for validating escape sequences in string literals.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class EscapeUtils {

    private static final String ESCAPABLE_CHARS = "btnfrs\"'\\$\r\n";

    private static final String SLASH_STR = "/";

    /**
     * Get the list of errors for each invalid escape sequence in a
     * string literal.
     *
     * The lexer accepts invalid escapes (see `InvalidEscape` in the grammar)
     * so that the error can be reported here at the exact position, instead
     * of causing the enclosing string to fail and be reported at the start
     * of the string.
     *
     * @param text
     * @param ctx
     */
    public static List<SyntaxException> checkEscapes(String text, ParserRuleContext ctx) {
        if( text.indexOf('\\') == -1 || text.startsWith(SLASH_STR) )
            return Collections.emptyList();
        var errors = new ArrayList<SyntaxException>();
        var line = ctx.getStart().getLine();
        var column = ctx.getStart().getCharPositionInLine() + 1;
        for( int i = 0; i < text.length(); i++ ) {
            var c = text.charAt(i);
            if( c == '\n' ) {
                line++;
                column = 1;
            }
            else if( c == '\\' && i + 1 < text.length() ) {
                var next = text.charAt(i + 1);
                if( !isEscapable(text, i + 1) )
                    errors.add(new SyntaxException("Invalid escape sequence: '\\" + next + "'", line, column, line, column + 2));
                i++;
                if( next == '\n' || next == '\r' ) {
                    line++;
                    column = 1;
                }
                else {
                    column += 2;
                }
            }
            else {
                column++;
            }
        }
        return errors;
    }

    private static boolean isEscapable(String text, int i) {
        var c = text.charAt(i);
        if( c == 'u' )
            return isUnicodeEscape(text, i + 1);
        return ESCAPABLE_CHARS.indexOf(c) != -1 || ('0' <= c && c <= '7');
    }

    /**
     * Determine whether a unicode escape is followed by four hex digits.
     * Groovy allows one or more `u`s after the backslash.
     */
    private static boolean isUnicodeEscape(String text, int i) {
        while( i < text.length() && text.charAt(i) == 'u' )
            i++;
        if( i + 4 > text.length() )
            return false;
        for( int j = i; j < i + 4; j++ ) {
            if( Character.digit(text.charAt(j), 16) == -1 )
                return false;
        }
        return true;
    }

}
