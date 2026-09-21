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
package nextflow.script.formatter;

import org.codehaus.groovy.ast.ASTNode;

/**
 * A comment in the source, collected from the token stream.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class Comment {

    /**
     * The comment text, with line endings normalized to LF.
     */
    public final String text;

    public final int line;
    public final int column;
    public final int lastLine;
    public final int lastColumn;

    /**
     * Whether the comment is the first thing on its line, as opposed
     * to an end-of-line comment after some code.
     */
    public final boolean ownLine;

    public Comment(String text, int line, int column, boolean ownLine) {
        this.text = text;
        this.line = line;
        this.column = column;
        this.ownLine = ownLine;

        var lastIndex = text.lastIndexOf('\n');
        this.lastLine = line + countNewLines(text);
        this.lastColumn = lastIndex < 0
            ? column + text.length()
            : text.length() - lastIndex;
    }

    public long start() {
        return position(line, column);
    }

    public long end() {
        return position(lastLine, lastColumn);
    }

    @Override
    public String toString() {
        return text + " (line " + line + ")";
    }

    /**
     * Pack a (line, column) pair into a single comparable value.
     */
    public static long position(int line, int column) {
        return ((long) line << 32) | (column & 0xffffffffL);
    }

    public static long start(ASTNode node) {
        return position(node.getLineNumber(), node.getColumnNumber());
    }

    public static long end(ASTNode node) {
        return position(node.getLastLineNumber(), node.getLastColumnNumber());
    }

    public static int countNewLines(String text) {
        int result = 0;
        for( int i = 0; i < text.length(); i++ ) {
            if( text.charAt(i) == '\n' )
                result++;
        }
        return result;
    }

    public static boolean hasPosition(ASTNode node) {
        return node != null && node.getLineNumber() > 0;
    }

    /**
     * Normalize a comment to a form that is invariant under re-indentation
     * and line re-wrapping, so that the comments of two sources can be
     * compared for equality.
     */
    public static String normalize(String text) {
        return text.replaceAll("\\s+", " ").strip();
    }

}
