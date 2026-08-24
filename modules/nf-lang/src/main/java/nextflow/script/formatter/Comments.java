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

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import nextflow.config.parser.ConfigLexer;
import nextflow.script.parser.ScriptLexer;
import org.antlr.v4.runtime.CharStreams;
import org.antlr.v4.runtime.CommonTokenStream;
import org.antlr.v4.runtime.Lexer;
import org.antlr.v4.runtime.Token;

/**
 * The comments in a source file, collected from the token stream
 * produced by the parse.
 *
 * Comments are lexed as newline tokens (see {@code ML_COMMENT},
 * {@code SL_COMMENT} and {@code SH_COMMENT} in the lexer grammars),
 * so they are identified by their leading characters -- a newline
 * token is either a line terminator or a comment.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class Comments {

    public static final Comments EMPTY = new Comments(Collections.emptyList(), new BitSet(), Collections.emptyMap(), null);

    private final List<Comment> comments;

    /**
     * The lines that contain at least one token, used to derive the
     * number of blank lines above a given line.
     */
    private final BitSet contentLines;

    /**
     * The column of the first token on each line, used to tell whether a
     * node begins its line -- i.e. whether its column is an indentation.
     */
    private final Map<Integer,Integer> firstColumns;

    private final String shebang;

    private Comments(List<Comment> comments, BitSet contentLines, Map<Integer,Integer> firstColumns, String shebang) {
        this.comments = comments;
        this.contentLines = contentLines;
        this.firstColumns = firstColumns;
        this.shebang = shebang;
    }

    /**
     * Determine whether a position is the first token on its line.
     *
     * @param line   1-based
     * @param column 1-based
     */
    public boolean startsLine(int line, int column) {
        var first = firstColumns.get(line);
        return first != null && first == column;
    }

    public List<Comment> getComments() {
        return comments;
    }

    public String getShebang() {
        return shebang;
    }

    public int blankLinesBefore(int line) {
        int result = 0;
        for( int i = line - 1; i >= 1 && !contentLines.get(i); i-- )
            result++;
        return result;
    }

    /**
     * Collect the comments from the token buffer built by the parse.
     *
     * @param tokens
     */
    public static Comments collect(List<Token> tokens) {
        var comments = new ArrayList<Comment>();
        var contentLines = new BitSet();
        var firstColumns = new HashMap<Integer,Integer>();
        String shebang = null;
        boolean first = true;
        int lastCodeLine = -1;

        for( var token : tokens ) {
            if( token.getType() == Token.EOF )
                continue;
            var text = token.getText();
            if( text == null || text.isEmpty() )
                continue;
            if( text.charAt(0) == '\n' || text.charAt(0) == '\r' )
                continue;

            text = text.replace("\r\n", "\n").replace('\r', '\n').stripTrailing();
            var line = token.getLine();
            var lastLine = line + Comment.countNewLines(text);
            contentLines.set(line, lastLine + 1);
            // the tokens are in source order, so the first token seen on a
            // line is the leftmost one
            firstColumns.putIfAbsent(line, token.getCharPositionInLine() + 1);

            if( isComment(text) ) {
                if( first && text.startsWith("#!") ) {
                    // the shebang must stay on the first line, so it is not
                    // subject to comment attachment
                    shebang = text;
                }
                else {
                    var column = token.getCharPositionInLine() + 1;
                    comments.add(new Comment(text, line, column, lastCodeLine != line));
                }
            }
            else {
                lastCodeLine = lastLine;
            }
            first = false;
        }

        return new Comments(comments, contentLines, firstColumns, shebang);
    }

    /**
     * Get the normalized text of every comment in a source string, sorted, so
     * that the comments of two sources can be compared for equality. Used to
     * verify that formatting did not lose or alter a comment.
     *
     * @param source
     * @param config whether the source is a config file rather than a script
     */
    public static List<String> textsOf(String source, boolean config) {
        Lexer lexer = config
            ? new ConfigLexer(CharStreams.fromString(source))
            : new ScriptLexer(CharStreams.fromString(source));
        lexer.removeErrorListeners();
        var stream = new CommonTokenStream(lexer);
        stream.fill();
        var comments = collect(stream.getTokens());
        var result = new ArrayList<String>();
        if( comments.getShebang() != null )
            result.add(Comment.normalize(comments.getShebang()));
        for( var comment : comments.getComments() )
            result.add(Comment.normalize(comment.text));
        Collections.sort(result);
        return result;
    }

    private static boolean isComment(String text) {
        return text.startsWith("//") || text.startsWith("/*") || text.startsWith("#!");
    }

}
