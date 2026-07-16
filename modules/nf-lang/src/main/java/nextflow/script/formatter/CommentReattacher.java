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
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import nextflow.config.ast.ConfigApplyBlockNode;
import nextflow.config.ast.ConfigAssignNode;
import nextflow.config.ast.ConfigBlockNode;
import nextflow.config.ast.ConfigIncludeNode;
import nextflow.config.ast.ConfigNode;
import nextflow.config.parser.ConfigLexer;
import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.ast.FunctionNode;
import nextflow.script.ast.OutputBlockNode;
import nextflow.script.ast.OutputNode;
import nextflow.script.ast.ParamBlockNode;
import nextflow.script.ast.ProcessNodeV1;
import nextflow.script.ast.ProcessNodeV2;
import nextflow.script.ast.RecordNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.WorkflowNode;
import nextflow.script.parser.ScriptLexer;
import org.antlr.v4.runtime.CharStreams;
import org.antlr.v4.runtime.CommonTokenStream;
import org.antlr.v4.runtime.Lexer;
import org.antlr.v4.runtime.Token;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.CodeVisitorSupport;
import org.codehaus.groovy.ast.ModuleNode;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.EmptyExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.expr.BinaryExpression;
import org.codehaus.groovy.ast.expr.GStringExpression;
import org.codehaus.groovy.ast.expr.ListExpression;
import org.codehaus.groovy.ast.expr.MapExpression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.NamedArgumentListExpression;
import org.codehaus.groovy.ast.expr.TupleExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.EmptyStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.IfStatement;
import org.codehaus.groovy.ast.stmt.ReturnStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.ast.stmt.TryCatchStatement;

/**
 * Re-derive comment metadata for an AST from the original source text,
 * so that the formatter can preserve every comment.
 *
 * The parser ({@code ScriptAstBuilder} / {@code ConfigAstBuilder})
 * only attaches comments that immediately precede a statement or declaration
 * (as {@link ASTNodeMarker#LEADING_COMMENTS}) and drops everything else:
 * comments after the last statement of a block, comments at the end of the
 * file, comments in empty blocks, comments between process/workflow sections,
 * and so on. It also stores trailing {@code \r\n} sequences verbatim, which
 * corrupts formatter output for files with Windows line endings.
 *
 * This class re-lexes the source (using the same ANTLR lexer as the parser,
 * so string literals and slashy strings are handled exactly), assigns every
 * comment token to an AST node based on its position, and rewrites the
 * comment metadata:
 *
 * <ul>
 *   <li>{@link ASTNodeMarker#LEADING_COMMENTS}: comments (and blank lines)
 *       preceding a node, in reverse source order (parser convention);</li>
 *   <li>{@link ASTNodeMarker#TRAILING_COMMENT}: a comment on the same line,
 *       after the end of a node;</li>
 *   <li>{@link #DANGLING_COMMENTS}: comments at the end of a block or file,
 *       after the last statement, in source order — emitted before the
 *       closing brace (or at end of file);</li>
 *   <li>{@link #DANGLING_AFTER}: comments after a node at the end of a
 *       section that has no block emission point (e.g. the last workflow
 *       take), in source order — emitted immediately after the node's
 *       line.</li>
 * </ul>
 *
 * All entries are normalized to LF line endings.
 *
 * @see Formatter
 * @see ScriptFormattingVisitor
 * @see nextflow.config.formatter.ConfigFormattingVisitor
 */
public class CommentReattacher {

    /**
     * Metadata key for comments at the end of a block-like region (list of
     * String in source order; blank lines are "\n" entries). Emitted before
     * the closing brace of the block, or at the end of the file.
     */
    public static final Object DANGLING_COMMENTS = Marker.DANGLING_COMMENTS;

    /**
     * Metadata key for comments that follow a node at the end of a section
     * that has no block emission point (list of String in source order).
     * Emitted immediately after the node's line.
     */
    public static final Object DANGLING_AFTER = Marker.DANGLING_AFTER;

    /**
     * Metadata key for a source region that is excluded from formatting
     * with `// fmt: skip` or `// fmt: off` ... `// fmt: on` (the verbatim
     * source text, emitted in place of the node).
     */
    public static final Object VERBATIM_SOURCE = Marker.VERBATIM_SOURCE;

    /**
     * Metadata key marking a node that is part of a verbatim region emitted
     * by an earlier node (the node itself emits nothing).
     */
    public static final Object VERBATIM_SUPPRESSED = Marker.VERBATIM_SUPPRESSED;

    private enum Marker {
        DANGLING_COMMENTS,
        DANGLING_AFTER,
        VERBATIM_SOURCE,
        VERBATIM_SUPPRESSED
    }

    /**
     * Rewrite the comment metadata of a module from its source text.
     *
     * @param moduleNode the AST (a {@link ScriptNode} or {@link ConfigNode})
     * @param sourceText the source text the AST was parsed from
     */
    public static void apply(ModuleNode moduleNode, String sourceText) {
        if( moduleNode instanceof ScriptNode sn ) {
            var source = lex(sourceText, false);
            var builder = new ScriptRegionBuilder();
            var root = builder.build(sn);
            new CommentReattacher(source, sourceText, builder.sectionHeads).run(root);
        }
        else if( moduleNode instanceof ConfigNode cn ) {
            var source = lex(sourceText, true);
            var root = new ConfigRegionBuilder().build(cn);
            new CommentReattacher(source, sourceText, Set.of()).run(root);
        }
    }

    /**
     * Count the comments in a source text. Used as a safety check that
     * formatting neither removes nor duplicates comments. The shebang line
     * is not counted (it is stored separately from comments).
     *
     * @param sourceText
     * @param configFile whether to lex as a config file instead of a script
     */
    public static int countComments(String sourceText, boolean configFile) {
        int count = 0;
        for( var token : lex(sourceText, configFile).tokens() ) {
            if( token.comment() && !token.shebang() )
                count++;
        }
        return count;
    }

    /// LEXING

    private record NlToken(
        String text,
        long start,         // encoded position of first character
        long end,           // encoded position just after last character
        int line,           // 1-based first line
        int lastLine,       // 1-based last line
        boolean comment     // false for pure newline tokens
    ) {
        boolean shebang() {
            return line == 1 && text.startsWith("#!");
        }

        boolean multiLine() {
            return line != lastLine;
        }

        /**
         * The comment text, normalized to LF line endings.
         */
        String normalized() {
            return text.replace("\r\n", "\n").replace("\r", "\n");
        }
    }

    /**
     * The comment and newline tokens of a source text, together with the
     * set of lines that contain code (i.e. any token other than a comment
     * or newline) and the set of lines that contain comments.
     */
    private record TokenizedSource(
        List<NlToken> tokens,
        Set<Integer> contentLines,
        Set<Integer> commentLines
    ) {}

    /**
     * Extract all comments and newlines from the source text. In the
     * script/config grammars, comments and newlines are both emitted as
     * {@code NL} tokens on the default channel, so lexing the source with
     * the real lexer recovers every comment with exact positions and no
     * false positives from string literals.
     */
    private static TokenizedSource lex(String sourceText, boolean configFile) {
        Lexer lexer = configFile
            ? new ConfigLexer(CharStreams.fromString(sourceText))
            : new ScriptLexer(CharStreams.fromString(sourceText));
        lexer.removeErrorListeners();
        int nlType = configFile ? ConfigLexer.NL : ScriptLexer.NL;

        var stream = new CommonTokenStream(lexer);
        stream.fill();

        var tokens = new ArrayList<NlToken>();
        var contentLines = new HashSet<Integer>();
        var commentLines = new HashSet<Integer>();
        for( var token : stream.getTokens() ) {
            var text = token.getText();
            var line = token.getLine();
            int newlines = 0;
            int lastLineLength = 0;
            for( int i = 0; i < text.length(); i++ ) {
                if( text.charAt(i) == '\n' ) {
                    newlines++;
                    lastLineLength = 0;
                }
                else {
                    lastLineLength++;
                }
            }
            var lastLine = line + newlines;
            if( token.getType() != nlType || token.getChannel() != Token.DEFAULT_CHANNEL ) {
                if( token.getType() != Token.EOF ) {
                    for( int l = line; l <= lastLine; l++ )
                        contentLines.add(l);
                }
                continue;
            }
            var col = token.getCharPositionInLine() + 1;
            var endCol = newlines > 0 ? lastLineLength + 1 : col + text.length();
            var comment = !text.isBlank();
            if( comment ) {
                for( int l = line; l <= lastLine; l++ )
                    commentLines.add(l);
            }
            tokens.add(new NlToken(text, pos(line, col), pos(lastLine, endCol), line, lastLine, comment));
        }
        return new TokenizedSource(tokens, contentLines, commentLines);
    }

    /// REGION MODEL

    /**
     * A region is a range of source text with a list of "anchor" nodes
     * (statements, parameters, declarations, ...) that the formatter emits
     * line-by-line, and a list of nested child regions (blocks, sections,
     * closures). Comments are assigned to the anchors and regions based on
     * their position.
     *
     * A container node (workflow, process, config block, ...) appears twice:
     * as an anchor of its parent region (so it can receive leading and
     * trailing comments) and as the holder of its own interior region (so
     * it can receive dangling comments, emitted before its closing brace).
     */
    /**
     * Which expression slots the anchors of a region provide for comments
     * inside multi-line expressions:
     *
     * NONE       -- none (e.g. directives, which are never wrapped);
     * ELEMENTS   -- elements of wrapped calls, lists and maps;
     * STATEMENTS -- elements and method-chain links (regular statements).
     */
    private enum SlotMode {
        NONE,
        ELEMENTS,
        STATEMENTS
    }

    private static class Region {
        ASTNode holder;         // receives DANGLING_COMMENTS
        long start;
        long end;
        SlotMode slotMode = SlotMode.NONE;
        List<ASTNode> anchors = new ArrayList<>();
        List<Region> children = new ArrayList<>();

        Region(ASTNode holder, long start, long end) {
            this.holder = holder;
            this.start = start;
            this.end = end;
        }

        Region(ASTNode holder, long start, long end, SlotMode slotMode) {
            this(holder, start, end);
            this.slotMode = slotMode;
        }

        void addAnchor(ASTNode node) {
            if( isPositioned(node) )
                anchors.add(node);
        }

        void addChild(Region region) {
            children.add(region);
        }

        void sortRecursive() {
            anchors.sort(Comparator.comparingLong(node -> startOf(node)));
            children.sort(Comparator.comparingLong(region -> region.start));
            for( var child : children )
                child.sortRecursive();
        }
    }

    /// COMMENT ASSIGNMENT

    private final List<NlToken> tokens;

    private final Set<Integer> contentLines;

    private final Set<Integer> commentLines;

    /**
     * Anchors whose source position includes their section label (e.g. the
     * process `when:` expression). The formatter emits its own blank line
     * before such sections, so the blank prefix of their leading comments
     * is trimmed to keep formatting stable.
     */
    private final Set<ASTNode> sectionHeads;

    private final String[] sourceLines;

    private final boolean[] claimed;

    private CommentReattacher(TokenizedSource source, String sourceText, Set<ASTNode> sectionHeads) {
        this.tokens = source.tokens();
        this.contentLines = source.contentLines();
        this.commentLines = source.commentLines();
        this.sectionHeads = sectionHeads;
        this.sourceLines = sourceText.replace("\r\n", "\n").replace("\r", "\n").split("\n", -1);
        this.claimed = new boolean[tokens.size()];
        // claim the shebang and its line terminator -- the shebang is
        // stored separately from comments and emitted by the formatter
        // together with its own line break
        for( int i = 0; i < tokens.size(); i++ ) {
            if( tokens.get(i).shebang() ) {
                claimed[i] = true;
                if( i + 1 < tokens.size() && !tokens.get(i + 1).comment() && tokens.get(i + 1).line() == tokens.get(i).line() )
                    claimed[i + 1] = true;
            }
        }
    }

    private void run(Region root) {
        root.sortRecursive();
        clearMetadata(root);
        applyFormatterDirectives(root);
        assignRegion(root);
    }

    /// FORMATTER DIRECTIVES (fmt: skip / fmt: off / fmt: on)

    private static final java.util.regex.Pattern FMT_DIRECTIVE = java.util.regex.Pattern.compile("//\\s*fmt:\\s*(skip|off|on)\\s*");

    private static String fmtDirective(NlToken token) {
        if( !token.comment() || token.multiLine() )
            return null;
        var matcher = FMT_DIRECTIVE.matcher(token.normalized());
        return matcher.matches() ? matcher.group(1) : null;
    }

    /**
     * Apply `// fmt: skip` and `// fmt: off` ... `// fmt: on` directives:
     * the statements and declarations they cover are emitted verbatim from
     * the source text instead of being formatted.
     *
     * A `fmt: skip` comment applies to the statement that ends on its line.
     * A `fmt: off` comment applies from its own line until the line of the
     * next `fmt: on` comment (or the end of the file), and should enclose
     * whole statements or declarations.
     */
    private void applyFormatterDirectives(Region root) {
        // determine the line ranges excluded from formatting
        var ranges = new ArrayList<int[]>();
        int rangeStart = -1;
        for( var token : tokens ) {
            var directive = fmtDirective(token);
            if( directive == null )
                continue;
            if( rangeStart >= 0 ) {
                if( "on".equals(directive) ) {
                    ranges.add(new int[] { rangeStart, token.line() });
                    rangeStart = -1;
                }
            }
            else if( "off".equals(directive) ) {
                rangeStart = token.line();
            }
            else if( "skip".equals(directive) ) {
                ranges.add(new int[] { token.line(), token.line() });
            }
        }
        if( rangeStart >= 0 )
            ranges.add(new int[] { rangeStart, Integer.MAX_VALUE });
        if( ranges.isEmpty() )
            return;

        // collect all anchors in source order
        var anchors = new ArrayList<ASTNode>();
        collectAnchors(root, anchors);
        anchors.sort(Comparator.comparingLong(node -> startOf(node)));

        for( var range : ranges ) {
            // find the anchors that intersect the range, extending the
            // range to cover them completely
            var startLine = range[0];
            var endLine = range[1];
            var overlapping = new ArrayList<ASTNode>();
            for( var anchor : anchors ) {
                if( anchor.getLineNumber() <= endLine && anchor.getLastLineNumber() >= startLine )
                    overlapping.add(anchor);
            }
            // an anchor that merely encloses the range (e.g. the workflow
            // around a `fmt: skip` statement) is not part of it -- keep
            // enclosing anchors only when the range covers them entirely
            var intersecting = new ArrayList<ASTNode>();
            for( var anchor : overlapping ) {
                var coveredByRange = anchor.getLineNumber() >= range[0] && anchor.getLastLineNumber() <= range[1];
                if( !coveredByRange && containsAnother(anchor, overlapping) )
                    continue;
                intersecting.add(anchor);
            }
            if( intersecting.isEmpty() )
                continue;
            for( var anchor : intersecting ) {
                startLine = Math.min(startLine, anchor.getLineNumber());
                endLine = Math.max(endLine, anchor.getLastLineNumber());
            }
            endLine = Math.min(endLine, sourceLines.length);
            // drop trailing empty lines (e.g. past the end of the file)
            while( endLine > startLine && sourceLines[endLine - 1].isEmpty() )
                endLine--;

            // claim the comments and newlines within the range -- they are
            // part of the verbatim text
            for( int i = 0; i < tokens.size(); i++ ) {
                var token = tokens.get(i);
                if( token.line() >= startLine && token.line() <= endLine )
                    claimed[i] = true;
            }

            // attach the verbatim text to the first anchor and suppress
            // the others
            var text = String.join("\n", java.util.Arrays.copyOfRange(sourceLines, startLine - 1, endLine));
            intersecting.get(0).putNodeMetaData(VERBATIM_SOURCE, text);
            for( int i = 1; i < intersecting.size(); i++ )
                intersecting.get(i).putNodeMetaData(VERBATIM_SUPPRESSED, Boolean.TRUE);
        }
    }

    private void collectAnchors(Region region, List<ASTNode> anchors) {
        anchors.addAll(region.anchors);
        for( var child : region.children )
            collectAnchors(child, anchors);
    }

    private static boolean containsAnother(ASTNode anchor, List<ASTNode> anchors) {
        for( var other : anchors ) {
            if( other == anchor )
                continue;
            if( startOf(anchor) <= startOf(other) && endOf(other) <= endOf(anchor) )
                return true;
        }
        return false;
    }

    /**
     * Clear existing comment metadata so that the assignment below is the
     * single source of truth (and re-running it is idempotent).
     */
    private void clearMetadata(Region region) {
        clearNodeMetadata(region.holder);
        for( var anchor : region.anchors )
            clearNodeMetadata(anchor);
        for( var child : region.children )
            clearMetadata(child);
    }

    private void clearNodeMetadata(ASTNode node) {
        if( node == null )
            return;
        node.removeNodeMetaData(ASTNodeMarker.LEADING_COMMENTS);
        node.removeNodeMetaData(ASTNodeMarker.TRAILING_COMMENT);
        node.removeNodeMetaData(DANGLING_COMMENTS);
        node.removeNodeMetaData(DANGLING_AFTER);
    }

    /**
     * Assign the comments in a region to its anchors and child regions.
     * Child regions claim their tokens first; the remaining tokens are
     * sliced by the gaps between anchors and child regions.
     */
    private void assignRegion(Region region) {
        for( var child : region.children )
            assignRegion(child);

        // build the list of slots (anchors and child regions) in source order
        var slots = new ArrayList<Object>();
        {
            int i = 0;
            int j = 0;
            while( i < region.anchors.size() || j < region.children.size() ) {
                var anchorStart = i < region.anchors.size() ? startOf(region.anchors.get(i)) : Long.MAX_VALUE;
                var childStart = j < region.children.size() ? region.children.get(j).start : Long.MAX_VALUE;
                if( anchorStart <= childStart )
                    slots.add(region.anchors.get(i++));
                else
                    slots.add(region.children.get(j++));
            }
        }

        // assign the tokens in each gap between consecutive slots
        long gapStart = region.start;
        Object prev = null;
        Object trailingPrev = null;
        for( var slot : slots ) {
            var slotStart = slot instanceof Region r ? r.start : startOf((ASTNode) slot);
            var slotEnd = slot instanceof Region r ? r.end : endOf((ASTNode) slot);
            assignGap(region, prev, trailingPrev, slot, gapTokens(gapStart, slotStart));
            if( slot instanceof ASTNode anchor )
                assignInteriorComments(anchor, region.slotMode);
            gapStart = Math.max(gapStart, slotEnd);
            prev = slot;
            // an anchor stays the trailing-comment candidate across the
            // child regions nested within it (e.g. the branches of an
            // if/else statement)
            if( slot instanceof ASTNode anchor )
                trailingPrev = anchor;
            else if( slot instanceof Region r && !(trailingPrev instanceof ASTNode anchor && r.end <= endOf(anchor)) )
                trailingPrev = r;
        }
        assignGap(region, prev, trailingPrev, null, gapTokens(gapStart, region.end));
    }

    /**
     * Collect the unclaimed tokens in a position range, and mark them as
     * claimed.
     */
    private List<NlToken> gapTokens(long start, long end) {
        var result = new ArrayList<NlToken>();
        for( int i = 0; i < tokens.size(); i++ ) {
            var token = tokens.get(i);
            if( claimed[i] || token.shebang() )
                continue;
            if( token.start() >= start && token.start() < end ) {
                claimed[i] = true;
                result.add(token);
            }
        }
        return result;
    }

    /**
     * Assign the tokens in the gap between two slots:
     *
     * <ul>
     *   <li>a single-line comment on the same line as the end of the
     *       previous slot becomes its trailing comment;</li>
     *   <li>the remaining tokens are split at the last newline that
     *       terminates a line of code (e.g. a section label or a closing
     *       brace): comments before it hang off the end of the previous
     *       slot, comments after it lead the next slot.</li>
     * </ul>
     */
    private void assignGap(Region region, Object prev, Object trailingPrev, Object next, List<NlToken> gap) {
        if( gap.isEmpty() )
            return;

        // claim a single-line comment on the same line as the end of the
        // previous slot as its trailing comment
        var trailingTarget = trailingTargetOf(trailingPrev);
        if( trailingTarget != null ) {
            var first = gap.get(0);
            if( first.comment() && !first.multiLine()
                    && first.line() == trailingTarget.getLastLineNumber()
                    && trailingTarget.getNodeMetaData(ASTNodeMarker.TRAILING_COMMENT) == null ) {
                trailingTarget.putNodeMetaData(ASTNodeMarker.TRAILING_COMMENT, first.normalized());
                gap = gap.subList(1, gap.size());
            }
        }

        // the trailing gap of the region dangles as a whole
        if( next == null ) {
            if( containsComment(gap) )
                appendRegionDangling(region, gap);
            return;
        }

        // split the gap at the last newline that terminates a line of code
        // (e.g. a section label): tokens after it lead the next slot
        int boundary = -1;
        for( int i = 0; i < gap.size(); i++ ) {
            var token = gap.get(i);
            if( !token.comment() && contentLines.contains(token.line()) )
                boundary = i;
        }
        var pre = gap.subList(0, boundary + 1);
        var post = new ArrayList<NlToken>(gap.subList(boundary + 1, gap.size()));

        if( containsComment(pre) ) {
            // comments adjacent to the next slot (after the last blank line)
            // also lead the next slot; the rest hang off the previous slot
            int split = 0;
            for( int i = pre.size() - 1; i >= 0; i-- ) {
                if( isBlankLine(pre.get(i)) ) {
                    split = i + 1;
                    break;
                }
            }
            var head = pre.subList(0, split);
            var tail = pre.subList(split, pre.size());
            if( !containsComment(head) ) {
                // the head is just the blank separator above the next slot
                post.addAll(0, tail);
            }
            else if( prev == null ) {
                // no previous slot to hang off of
                post.addAll(0, pre);
            }
            else {
                if( prev instanceof Region prevRegion )
                    appendRegionDangling(prevRegion, head);
                else
                    appendDangling((ASTNode) prev, DANGLING_AFTER, head);
                post.addAll(0, tail);
            }
        }

        // strip blank lines at the start of a block or section (the
        // formatter provides canonical spacing there); blank lines between
        // two statements are preserved
        if( !(prev instanceof ASTNode && next instanceof ASTNode) ) {
            while( !post.isEmpty() && isBlankLine(post.get(0)) )
                post.remove(0);
        }

        // assign the remaining tokens (comments and blank lines) as the
        // leading comments of the next slot; a blank-only gap before a
        // section is just the section separator, which the formatter
        // re-creates itself
        if( next instanceof Region ) {
            if( !containsComment(post) )
                return;
            // comments separated from the section by a blank line hang off
            // the previous slot instead
            int split = 0;
            for( int i = post.size() - 1; i >= 0; i-- ) {
                if( isBlankLine(post.get(i)) ) {
                    split = i + 1;
                    break;
                }
            }
            var head = post.subList(0, split);
            if( prev != null && containsComment(head) ) {
                if( prev instanceof Region prevRegion )
                    appendRegionDangling(prevRegion, head);
                else
                    appendDangling((ASTNode) prev, DANGLING_AFTER, head);
                post = new ArrayList<NlToken>(post.subList(split, post.size()));
            }
            if( !containsComment(post) )
                return;
        }
        var target = next instanceof ASTNode anchor ? anchor : firstAnchorOf((Region) next);
        if( target != null ) {
            if( sectionHeads.contains(target) ) {
                while( !post.isEmpty() && !post.get(0).comment() )
                    post.remove(0);
            }
            setLeadingComments(target, post);
        }
        else if( containsComment(post) )
            appendRegionDangling((Region) next, post);
    }

    /**
     * Assign dangling comments to a region: to its holder if it has one,
     * otherwise after its last anchor.
     */
    private void appendRegionDangling(Region region, List<NlToken> gap) {
        if( region.holder != null )
            appendDangling(region.holder, DANGLING_COMMENTS, gap);
        else if( !region.anchors.isEmpty() )
            appendDangling(region.anchors.get(region.anchors.size() - 1), DANGLING_AFTER, gap);
        else if( !region.children.isEmpty() )
            appendRegionDangling(region.children.get(region.children.size() - 1), gap);
    }

    /**
     * The node that can claim a trailing comment after the given slot:
     * an anchor claims it directly; for a region, the holder claims it if
     * it is a container node (e.g. a comment after the closing brace of a
     * workflow), otherwise the last anchor of the region.
     */
    private ASTNode trailingTargetOf(Object slot) {
        if( slot instanceof ASTNode node )
            return node;
        if( slot instanceof Region region ) {
            if( region.holder != null && !(region.holder instanceof BlockStatement) && isPositioned(region.holder) )
                return region.holder;
            if( !region.anchors.isEmpty() ) {
                var lastAnchor = region.anchors.get(region.anchors.size() - 1);
                var lastChild = region.children.isEmpty() ? null : trailingTargetOf(region.children.get(region.children.size() - 1));
                if( lastChild != null && endOf(lastChild) > endOf(lastAnchor) )
                    return lastChild;
                return lastAnchor;
            }
            if( !region.children.isEmpty() )
                return trailingTargetOf(region.children.get(region.children.size() - 1));
        }
        return null;
    }

    /**
     * The first anchor of a region, descending into child regions.
     */
    private ASTNode firstAnchorOf(Region region) {
        var firstAnchor = region.anchors.isEmpty() ? null : region.anchors.get(0);
        var firstChild = region.children.isEmpty() ? null : firstAnchorOf(region.children.get(0));
        if( firstAnchor == null )
            return firstChild;
        if( firstChild == null )
            return firstAnchor;
        return startOf(firstChild) < startOf(firstAnchor) ? firstChild : firstAnchor;
    }

    /**
     * Convert gap tokens to metadata entries: comments keep their
     * (normalized) text, newlines on blank lines become "\n" entries, and
     * newlines that terminate a line of code are dropped (the formatter
     * re-creates those line breaks itself).
     */
    private List<String> gapEntries(List<NlToken> gap) {
        var entries = new ArrayList<String>();
        var prevBlank = false;
        for( var token : gap ) {
            if( token.comment() ) {
                entries.add(token.normalized());
                prevBlank = false;
                // a comment that shares its line with code has no newline
                // entry of its own (line terminators of code lines are
                // dropped) -- give it one so that whatever the formatter
                // emits next does not end up on the comment line
                if( contentLines.contains(token.lastLine()) )
                    entries.add("\n");
            }
            else if( contentLines.contains(token.line()) ) {
                // line terminator of a line of code -- the formatter
                // re-creates the line break itself
            }
            else if( commentLines.contains(token.line()) ) {
                // line terminator of a comment line
                entries.add("\n");
                prevBlank = false;
            }
            else {
                // blank line -- collapse consecutive blank lines into one
                if( !prevBlank )
                    entries.add("\n");
                prevBlank = true;
            }
        }
        return entries;
    }

    /**
     * Store the leading comments of a node, following the parser convention:
     * entries in reverse source order.
     */
    private void setLeadingComments(ASTNode node, List<NlToken> gap) {
        var entries = gapEntries(gap);
        if( entries.isEmpty() )
            return;
        var reversed = new ArrayList<String>(entries.size());
        for( int i = entries.size() - 1; i >= 0; i-- )
            reversed.add(entries.get(i));
        var existing = (List<String>) node.getNodeMetaData(ASTNodeMarker.LEADING_COMMENTS);
        if( existing != null ) {
            // entries stored earlier come later in the source, so the new
            // (earlier) entries are appended to the reversed list
            existing.addAll(reversed);
        }
        else {
            node.putNodeMetaData(ASTNodeMarker.LEADING_COMMENTS, reversed);
        }
    }

    /**
     * Append dangling comments to a node under the given key: entries in
     * source order, with blank lines trimmed from the end (blank lines
     * before the first comment are preserved).
     */
    private void appendDangling(ASTNode node, Object key, List<NlToken> gap) {
        var entries = gapEntries(gap);
        while( !entries.isEmpty() && "\n".equals(entries.get(entries.size() - 1)) )
            entries.remove(entries.size() - 1);
        if( entries.isEmpty() || !entries.stream().anyMatch(entry -> !"\n".equals(entry)) )
            return;

        var existing = (List<String>) node.getNodeMetaData(key);
        if( existing != null ) {
            existing.addAll(entries);
        }
        else {
            node.putNodeMetaData(key, entries);
        }
    }

    /**
     * Assign comments inside an anchor (e.g. inside a multi-line expression)
     * that were not claimed by a child region. Comments that fall on an
     * expression slot -- before a link of the root method chain, or before
     * an element of a wrapped call, list or map -- are attached to that
     * slot and emitted in place; the rest are hoisted into the anchor's
     * leading comments rather than being dropped (the formatter cannot emit
     * comments inside arbitrary expressions).
     */
    private void assignInteriorComments(ASTNode anchor, SlotMode mode) {
        var interior = gapTokens(startOf(anchor), endOf(anchor));
        if( !containsComment(interior) )
            return;

        var slots = collectExpressionSlots(anchor, mode);
        var attached = new java.util.LinkedHashMap<ASTNode, List<String>>();
        var leftover = new ArrayList<NlToken>();
        for( var token : interior ) {
            if( !token.comment() )
                continue;
            var slot = findSlot(slots, token);
            if( slot != null )
                attached.computeIfAbsent(slot.target(), k -> new ArrayList<>()).add(token.normalized());
            else
                leftover.add(token);
        }

        // store slot comments following the leading comments convention:
        // entries in reverse source order, each comment followed by its
        // line terminator
        for( var slotEntry : attached.entrySet() ) {
            var entries = new ArrayList<String>();
            for( var comment : slotEntry.getValue() ) {
                entries.add(0, "\n");
                entries.add(0, comment);
            }
            var reversed = new ArrayList<String>(entries.size());
            for( int i = entries.size() - 1; i >= 0; i-- )
                reversed.add(entries.get(i));
            slotEntry.getKey().putNodeMetaData(ASTNodeMarker.LEADING_COMMENTS, reversed);
        }

        if( leftover.isEmpty() )
            return;
        var existing = (List<String>) anchor.getNodeMetaData(ASTNodeMarker.LEADING_COMMENTS);
        var entries = existing != null ? existing : new ArrayList<String>();
        for( var token : leftover ) {
            // prepend in reverse source order
            entries.add(0, "\n");
            entries.add(1, token.normalized());
        }
        if( existing == null )
            anchor.putNodeMetaData(ASTNodeMarker.LEADING_COMMENTS, entries);
    }

    /// EXPRESSION SLOTS

    /**
     * A position range where comments can be emitted inside an expression:
     * comments in [start, end) are attached to the target node and emitted
     * before it.
     */
    private record ExprSlot(
        long start,
        long end,
        ASTNode target
    ) {}

    private static ExprSlot findSlot(List<ExprSlot> slots, NlToken token) {
        for( var slot : slots ) {
            if( token.start() >= slot.start() && token.start() < slot.end() )
                return slot;
        }
        return null;
    }

    /**
     * Collect the expression slots of an anchor: the links of the root
     * method chain (which the formatter can wrap onto one line per link)
     * and the elements of multi-line calls, lists and maps (which the
     * formatter wraps onto one line per element).
     */
    private List<ExprSlot> collectExpressionSlots(ASTNode anchor, SlotMode mode) {
        if( mode == SlotMode.NONE )
            return List.of();

        var slots = new ArrayList<ExprSlot>();

        // the links of the root method chain -- only a chain at the root of
        // a regular statement can be wrapped by the formatter
        if( mode == SlotMode.STATEMENTS ) {
            var expr = rootExpressionOf(anchor);
            while( expr instanceof MethodCallExpression mce && !mce.isImplicitThis() ) {
                var receiver = mce.getObjectExpression();
                if( isPositioned(receiver) && isPositioned(mce.getMethod()) )
                    slots.add(new ExprSlot(endOf(receiver), startOf(mce.getMethod()), mce));
                expr = receiver;
            }
        }

        // the elements of multi-line calls, lists and maps
        // NOTE: config statements do not implement the Groovy code visitor,
        // so they must be matched before the Statement branch
        var collector = new ElementSlotCollector(slots);
        if( anchor instanceof ConfigAssignNode can )
            can.value.visit(collector);
        else if( anchor instanceof ConfigIncludeNode cin )
            cin.source.visit(collector);
        else if( anchor instanceof nextflow.script.ast.ParamNodeV1 pn )
            pn.value.visit(collector);
        else if( anchor instanceof Statement stmt )
            stmt.visit(collector);
        else if( anchor instanceof Parameter p && p.hasInitialExpression() )
            p.getInitialExpression().visit(collector);
        else if( anchor instanceof Expression e )
            e.visit(collector);

        slots.sort(Comparator.comparingLong(ExprSlot::start));
        return slots;
    }

    private static Expression rootExpressionOf(ASTNode anchor) {
        Expression expr = null;
        if( anchor instanceof ExpressionStatement es )
            expr = es.getExpression();
        else if( anchor instanceof ReturnStatement rs )
            expr = rs.getExpression();
        if( expr instanceof BinaryExpression be && be.getOperation().isA(org.codehaus.groovy.syntax.Types.ASSIGNMENT_OPERATOR) )
            expr = be.getRightExpression();
        return expr;
    }

    /**
     * Collect element slots for the multi-line calls, lists and maps in an
     * expression tree. Closure bodies are separate regions and GString
     * interpolations cannot hold comments, so they are not descended into.
     */
    private static class ElementSlotCollector extends CodeVisitorSupport {

        private final List<ExprSlot> slots;

        ElementSlotCollector(List<ExprSlot> slots) {
            this.slots = slots;
        }

        @Override
        public void visitMethodCallExpression(MethodCallExpression node) {
            if( node.getArguments() instanceof TupleExpression args && isMultiLine(args) ) {
                var elements = new ArrayList<Expression>();
                for( var arg : args.getExpressions() ) {
                    if( arg instanceof NamedArgumentListExpression named )
                        elements.addAll(named.getMapEntryExpressions());
                    else if( arg instanceof ClosureExpression )
                        continue;   // trailing closures are emitted outside the parentheses
                    else
                        elements.add(arg);
                }
                addElementSlots(endOf(node.getMethod()), elements);
            }
            super.visitMethodCallExpression(node);
        }

        @Override
        public void visitListExpression(ListExpression node) {
            if( isMultiLine(node) )
                addElementSlots(startOf(node), node.getExpressions());
            super.visitListExpression(node);
        }

        @Override
        public void visitMapExpression(MapExpression node) {
            if( !(node instanceof NamedArgumentListExpression) && isMultiLine(node) )
                addElementSlots(startOf(node), node.getMapEntryExpressions());
            super.visitMapExpression(node);
        }

        @Override
        public void visitClosureExpression(ClosureExpression node) {
            // do not descend -- closure bodies are separate regions
        }

        @Override
        public void visitGStringExpression(GStringExpression node) {
            // do not descend -- interpolations are emitted inline
        }

        @Override
        public void visitBlockStatement(BlockStatement node) {
            // do not descend -- nested blocks are separate regions
        }

        private void addElementSlots(long constructStart, List<? extends Expression> elements) {
            var sorted = new ArrayList<Expression>();
            for( var element : elements ) {
                if( isPositioned(element) )
                    sorted.add(element);
            }
            sorted.sort(Comparator.comparingLong(node -> startOf(node)));
            var gapStart = constructStart;
            for( var element : sorted ) {
                slots.add(new ExprSlot(gapStart, startOf(element), element));
                gapStart = endOf(element);
            }
        }

        private static boolean isMultiLine(Expression node) {
            return isPositioned(node) && node.getLineNumber() < node.getLastLineNumber();
        }
    }

    private static boolean containsComment(List<NlToken> gap) {
        for( var token : gap ) {
            if( token.comment() )
                return true;
        }
        return false;
    }

    /**
     * Determine whether a newline token terminates a blank line (a line
     * with no code and no comment).
     */
    private boolean isBlankLine(NlToken token) {
        return !token.comment()
            && !contentLines.contains(token.line())
            && !commentLines.contains(token.line());
    }

    /// POSITIONS

    private static long pos(int line, int col) {
        return (long) line * 1_000_000L + col;
    }

    private static long startOf(ASTNode node) {
        return pos(node.getLineNumber(), node.getColumnNumber());
    }

    private static long endOf(ASTNode node) {
        return pos(node.getLastLineNumber(), node.getLastColumnNumber());
    }

    private static boolean isPositioned(ASTNode node) {
        return node != null
            && !(node instanceof EmptyStatement)
            && !(node instanceof EmptyExpression)
            && node.getLineNumber() > 0
            && node.getLastLineNumber() > 0;
    }

    private static boolean isPositionedBlock(Statement node) {
        return node instanceof BlockStatement && isPositioned(node);
    }

    /// SCRIPT REGIONS

    private static class ScriptRegionBuilder {

        final Set<ASTNode> sectionHeads = new HashSet<>();

        Region build(ScriptNode scriptNode) {
            var root = new Region(scriptNode, 0, Long.MAX_VALUE, SlotMode.ELEMENTS);
            for( var decl : scriptNode.getDeclarations() ) {
                if( !isPositioned(decl) )
                    continue;
                root.addAnchor(decl);
                if( decl instanceof WorkflowNode wn )
                    root.addChild(workflowRegion(wn));
                else if( decl instanceof ProcessNodeV1 pn )
                    root.addChild(processRegionV1(pn));
                else if( decl instanceof ProcessNodeV2 pn )
                    root.addChild(processRegionV2(pn));
                else if( decl instanceof FunctionNode fn )
                    root.addChild(functionRegion(fn));
                else if( decl instanceof ParamBlockNode pbn )
                    root.addChild(paramsRegion(pbn));
                else if( decl instanceof OutputBlockNode obn )
                    root.addChild(outputsRegion(obn));
                else if( decl instanceof RecordNode rn )
                    root.addChild(recordRegion(rn));
                else if( decl instanceof ClassNode cn && cn.isEnum() )
                    root.addChild(enumRegion(cn));
            }
            return root;
        }

        private Region workflowRegion(WorkflowNode node) {
            var region = new Region(node, startOf(node), endOf(node), SlotMode.ELEMENTS);
            for( var take : node.getParameters() )
                region.addAnchor(take);
            addBlockRegion(region, node.main, SlotMode.STATEMENTS);
            addBlockRegion(region, node.emits, SlotMode.ELEMENTS);
            addBlockRegion(region, node.publishers, SlotMode.ELEMENTS);
            addBlockRegion(region, node.onComplete, SlotMode.STATEMENTS);
            addBlockRegion(region, node.onError, SlotMode.STATEMENTS);
            return region;
        }

        private Region processRegionV1(ProcessNodeV1 node) {
            var region = new Region(node, startOf(node), endOf(node), SlotMode.ELEMENTS);
            addBlockRegion(region, node.directives, SlotMode.NONE);
            addBlockRegion(region, node.inputs, SlotMode.NONE);
            addBlockRegion(region, node.outputs, SlotMode.NONE);
            if( isPositioned(node.when) ) {
                region.addAnchor(node.when);
                sectionHeads.add(node.when);
            }
            addBlockRegion(region, node.exec, SlotMode.STATEMENTS);
            addBlockRegion(region, node.stub, SlotMode.STATEMENTS);
            return region;
        }

        private Region processRegionV2(ProcessNodeV2 node) {
            var region = new Region(node, startOf(node), endOf(node), SlotMode.ELEMENTS);
            addBlockRegion(region, node.directives, SlotMode.NONE);
            for( var input : node.inputs )
                region.addAnchor(input);
            addBlockRegion(region, node.stagers, SlotMode.NONE);
            addBlockRegion(region, node.outputs, SlotMode.ELEMENTS);
            addBlockRegion(region, node.topics, SlotMode.STATEMENTS);
            if( isPositioned(node.when) ) {
                region.addAnchor(node.when);
                sectionHeads.add(node.when);
            }
            addBlockRegion(region, node.exec, SlotMode.STATEMENTS);
            addBlockRegion(region, node.stub, SlotMode.STATEMENTS);
            return region;
        }

        private Region functionRegion(FunctionNode node) {
            var region = new Region(node, startOf(node), endOf(node));
            addBlockRegion(region, node.getCode(), SlotMode.STATEMENTS);
            return region;
        }

        private Region paramsRegion(ParamBlockNode node) {
            var region = new Region(node, startOf(node), endOf(node), SlotMode.ELEMENTS);
            for( var param : node.declarations )
                region.addAnchor(param);
            return region;
        }

        private Region outputsRegion(OutputBlockNode node) {
            var region = new Region(node, startOf(node), endOf(node));
            for( var output : node.declarations ) {
                region.addAnchor(output);
                var outputRegion = outputRegion(output);
                if( outputRegion != null )
                    region.addChild(outputRegion);
            }
            return region;
        }

        private Region outputRegion(OutputNode node) {
            // an output declaration may not have a source position of its
            // own -- fall back to the position of its body
            var body = node.body instanceof BlockStatement block && isPositioned(block)
                ? (BlockStatement) node.body
                : null;
            var start = isPositioned(node) ? startOf(node) : body != null ? startOf(body) : -1;
            var end = isPositioned(node) ? endOf(node) : body != null ? endOf(body) : -1;
            if( start < 0 )
                return null;
            var region = new Region(node, start, end);
            if( body != null )
                region.addChild(blockRegion(body, startOf(body), endOf(body), SlotMode.NONE));
            return region;
        }

        private Region recordRegion(RecordNode node) {
            var region = new Region(node, startOf(node), endOf(node));
            for( var field : node.getFields() )
                region.addAnchor(field);
            return region;
        }

        private Region enumRegion(ClassNode node) {
            var region = new Region(node, startOf(node), endOf(node));
            for( var field : node.getFields() )
                region.addAnchor(field);
            return region;
        }

        /**
         * Add a child region for a block statement (a workflow/process
         * section, function body, ...) with its statements as anchors,
         * recursing into nested blocks and closures.
         */
        private void addBlockRegion(Region parent, Statement statement, SlotMode slotMode) {
            if( !isPositionedBlock(statement) )
                return;
            var block = (BlockStatement) statement;
            if( block.getStatements().isEmpty() )
                return;
            parent.addChild(blockRegion(block, startOf(block), endOf(block), slotMode));
        }

        private Region blockRegion(BlockStatement block, long start, long end) {
            return blockRegion(block, start, end, SlotMode.STATEMENTS);
        }

        private Region blockRegion(BlockStatement block, long start, long end, SlotMode slotMode) {
            var region = new Region(block, start, end, slotMode);
            for( var stmt : block.getStatements() ) {
                if( !isPositioned(stmt) )
                    continue;
                region.addAnchor(stmt);
                addStatementRegions(region, stmt);
            }
            return region;
        }

        /**
         * Add child regions for the nested blocks of a statement: if/else
         * branches, try/catch blocks, and closure bodies.
         */
        private void addStatementRegions(Region parent, Statement stmt) {
            if( stmt instanceof IfStatement is ) {
                addIfRegions(parent, is);
                return;
            }
            if( stmt instanceof TryCatchStatement tcs ) {
                if( isPositionedBlock(tcs.getTryStatement()) ) {
                    var tryBlock = (BlockStatement) tcs.getTryStatement();
                    var end = tcs.getCatchStatements().isEmpty()
                        ? endOf(tcs)
                        : startOf(tcs.getCatchStatements().get(0));
                    parent.addChild(blockRegion(tryBlock, startOf(tcs), end));
                }
                var catches = tcs.getCatchStatements();
                for( int i = 0; i < catches.size(); i++ ) {
                    var cs = catches.get(i);
                    if( !isPositionedBlock(cs.getCode()) )
                        continue;
                    var end = i + 1 < catches.size() ? startOf(catches.get(i + 1)) : endOf(tcs);
                    parent.addChild(blockRegion((BlockStatement) cs.getCode(), startOf(cs), end));
                }
                return;
            }
            // collect closures in the statement's expressions
            var collector = new ClosureCollector();
            stmt.visit(collector);
            for( var closure : collector.closures ) {
                if( closure.getCode() instanceof BlockStatement block && isPositioned(closure) )
                    parent.addChild(blockRegion(block, startOf(closure), endOf(closure)));
            }
        }

        private void addIfRegions(Region parent, IfStatement is) {
            var ifBlock = is.getIfBlock();
            var elseStmt = is.getElseBlock();
            var hasRealElse = elseStmt != null && !(elseStmt instanceof EmptyStatement);

            // an empty else block has no source position; comments after
            // the last statement of the if block are assigned to it
            long elseBoundary = endOf(is);
            if( hasRealElse ) {
                if( isPositioned(elseStmt) )
                    elseBoundary = startOf(elseStmt);
                else if( isPositionedBlock(ifBlock) )
                    elseBoundary = endOf(ifBlock);
                else
                    elseBoundary = startOf(is);
            }

            if( isPositionedBlock(ifBlock) )
                parent.addChild(blockRegion((BlockStatement) ifBlock, startOf(is), elseBoundary));

            if( elseStmt instanceof IfStatement chained ) {
                addIfRegions(parent, chained);
            }
            else if( hasRealElse && elseStmt instanceof BlockStatement elseBlock ) {
                var start = isPositioned(elseBlock) ? startOf(elseBlock) : elseBoundary;
                parent.addChild(blockRegion(elseBlock, start, endOf(is)));
            }
        }
    }

    /**
     * Collect the top-level closures of a statement, without descending
     * into nested closures (they are handled when the parent closure's
     * block region is built) or nested blocks (they are handled as child
     * regions of the statement).
     */
    private static class ClosureCollector extends CodeVisitorSupport {
        List<ClosureExpression> closures = new ArrayList<>();

        @Override
        public void visitClosureExpression(ClosureExpression node) {
            closures.add(node);
            // do not descend
        }

        @Override
        public void visitBlockStatement(BlockStatement node) {
            // do not descend
        }

        @Override
        public void visitIfElse(IfStatement node) {
            node.getBooleanExpression().visit(this);
        }

        @Override
        public void visitTryCatchFinally(TryCatchStatement node) {
        }
    }

    /// CONFIG REGIONS

    private static class ConfigRegionBuilder {

        Region build(ConfigNode configNode) {
            var root = new Region(configNode, 0, Long.MAX_VALUE, SlotMode.ELEMENTS);
            for( var stmt : configNode.getConfigStatements() )
                addConfigStatement(root, stmt);
            return root;
        }

        private void addConfigStatement(Region parent, ASTNode stmt) {
            if( !isPositioned(stmt) )
                return;
            parent.addAnchor(stmt);
            if( stmt instanceof ConfigBlockNode block ) {
                var region = new Region(block, startOf(block), endOf(block), SlotMode.ELEMENTS);
                for( var child : block.statements )
                    addConfigStatement(region, child);
                parent.addChild(region);
            }
            else if( stmt instanceof ConfigApplyBlockNode block ) {
                var region = new Region(block, startOf(block), endOf(block), SlotMode.NONE);
                for( var child : block.statements )
                    addConfigStatement(region, child);
                parent.addChild(region);
            }
            else {
                var value =
                    stmt instanceof ConfigAssignNode can ? can.value :
                    stmt instanceof ConfigIncludeNode cin ? cin.source :
                    null;
                if( value != null )
                    addExpressionRegions(parent, value);
            }
        }

        private void addExpressionRegions(Region parent, Expression expression) {
            var collector = new ClosureCollector();
            expression.visit(collector);
            for( var closure : collector.closures ) {
                if( closure.getCode() instanceof BlockStatement block && isPositioned(closure) )
                    parent.addChild(configClosureRegion(block, closure));
            }
        }

        private Region configClosureRegion(BlockStatement block, ClosureExpression closure) {
            var region = new Region(block, startOf(closure), endOf(closure), SlotMode.STATEMENTS);
            for( var stmt : block.getStatements() ) {
                if( !isPositioned(stmt) )
                    continue;
                region.addAnchor(stmt);
                var collector = new ClosureCollector();
                stmt.visit(collector);
                for( var nested : collector.closures ) {
                    if( nested.getCode() instanceof BlockStatement nestedBlock && isPositioned(nested) )
                        region.addChild(configClosureRegion(nestedBlock, nested));
                }
            }
            return region;
        }
    }

}
