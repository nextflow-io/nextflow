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
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Supplier;
import java.util.regex.Pattern;

import nextflow.config.ast.ConfigApplyBlockNode;
import nextflow.config.ast.ConfigAssignNode;
import nextflow.config.ast.ConfigBlockNode;
import nextflow.config.ast.ConfigIncludeNode;
import nextflow.config.ast.ConfigNode;
import nextflow.config.ast.ConfigStatement;
import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.ast.AgentNode;
import nextflow.script.ast.FunctionNode;
import nextflow.script.ast.OutputBlockNode;
import nextflow.script.ast.ProcessNodeV1;
import nextflow.script.ast.ProcessNodeV2;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.WorkflowNode;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.CodeVisitorSupport;
import org.codehaus.groovy.ast.ModuleNode;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.stmt.AssertStatement;
import org.codehaus.groovy.ast.stmt.EmptyStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.IfStatement;
import org.codehaus.groovy.ast.stmt.ReturnStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.ast.stmt.ThrowStatement;
import org.codehaus.groovy.ast.stmt.TryCatchStatement;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.runtime.IOGroovyMethods;

import static nextflow.script.ast.ASTUtils.*;
import static nextflow.script.formatter.Comment.*;

/**
 * Resolve `// fmt: skip` and `// fmt: off` ... `// fmt: on` directives into
 * the source ranges they exclude from formatting, following the semantics of
 * Black (language-server issue #75).
 *
 * Like {@link CommentAttacher}, this is derived fresh for each formatting
 * pass rather than stored in node metadata, since the same AST may be
 * formatted many times.
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
public class FmtDirectives {

    private static final FmtDirectives EMPTY = new FmtDirectives();

    private static final Pattern FMT_DIRECTIVE = Pattern.compile("//\\s*fmt:\\s*(skip|off|on)\\s*");

    private final Map<ASTNode,String> verbatimSource = new IdentityHashMap<>();

    /**
     * The first source line of the verbatim text assigned to an anchor --
     * usually the anchor's own line, but for `fmt: off` it is the line of the
     * `off` comment, which may precede the anchor's own leading comments.
     * Used by {@link Formatter#appendVerbatim} to avoid re-printing (or
     * double-counting the blank lines above) a leading comment that is
     * already part of the verbatim text.
     */
    private final Map<ASTNode,Integer> verbatimStart = new IdentityHashMap<>();

    private final Set<ASTNode> suppressed = Collections.newSetFromMap(new IdentityHashMap<>());

    /**
     * The comments whose own line falls within some verbatim range -- they
     * are printed as part of the raw text rather than individually, so the
     * normal comment-printing paths (dangling, trailing, a non-anchor's
     * leading comments) must not print them again.
     *
     * This is tracked separately from {@code CommentAttacher.markEmitted},
     * which is also (re-)called by the normal printing paths themselves --
     * including, harmlessly, more than once per comment when
     * {@code Formatter.emitWrappable} retries a statement with wrapping
     * enabled. Reusing that as a "already printed, skip it" signal would
     * make the first (rolled-back) attempt suppress the comment on retry.
     */
    private final Set<Comment> verbatimComments = Collections.newSetFromMap(new IdentityHashMap<>());

    private FmtDirectives() {
    }

    /**
     * Derive the verbatim regions for a module. As with {@link CommentAttacher},
     * the result belongs to a single formatting pass.
     *
     * The raw source text is supplied lazily and read only when the file
     * actually contains a fmt directive, so the common case (no directives)
     * never pays for reading the whole file.
     *
     * @param moduleNode
     * @param sourceText supplies the raw source text; may return null if it
     *                   could not be read, in which case directives are inactive
     * @param comments
     */
    public static FmtDirectives of(ModuleNode moduleNode, Supplier<String> sourceText, CommentAttacher comments) {
        if( moduleNode == null || sourceText == null )
            return EMPTY;
        var commentsMeta = (Comments) moduleNode.getNodeMetaData(ASTNodeMarker.COMMENTS);
        if( commentsMeta == null || commentsMeta.getComments().isEmpty() )
            return EMPTY;
        var ranges = computeRanges(commentsMeta.getComments());
        if( ranges.isEmpty() )
            return EMPTY;
        var text = sourceText.get();
        if( text == null )
            return EMPTY;
        var result = new FmtDirectives();
        result.apply(moduleNode, text, ranges, commentsMeta, comments);
        return result;
    }

    /**
     * Read the raw source text of a compilation unit, or null if it cannot be
     * read (in which case fmt directives are simply inactive).
     *
     * @param sourceUnit
     */
    public static String readSourceText(SourceUnit sourceUnit) {
        try {
            var reader = sourceUnit.getSource().getReader();
            return reader != null ? IOGroovyMethods.getText(reader) : null;
        }
        catch( Exception e ) {
            return null;
        }
    }

    public String verbatimText(ASTNode node) {
        return verbatimSource.get(node);
    }

    /**
     * The first source line covered by an anchor's verbatim text (see
     * {@link #verbatimStart}), or -1 if the anchor has none.
     */
    public int verbatimStartLine(ASTNode node) {
        var line = verbatimStart.get(node);
        return line != null ? line : -1;
    }

    public boolean isSuppressed(ASTNode node) {
        return suppressed.contains(node);
    }

    /**
     * Determine whether a comment's own line falls within some verbatim
     * range, i.e. it is already printed as part of a verbatim text.
     */
    public boolean coversComment(Comment comment) {
        return verbatimComments.contains(comment);
    }

    /**
     * Whether any comment at all is covered by a verbatim region. Lets the
     * common (no-directive) case skip the per-comment {@link #coversComment}
     * filtering entirely.
     */
    public boolean coversAnyComment() {
        return !verbatimComments.isEmpty();
    }

    /// RANGE RESOLUTION

    private void apply(ModuleNode moduleNode, String sourceText, List<int[]> ranges, Comments commentsMeta, CommentAttacher comments) {
        var sourceLines = sourceText.replace("\r\n", "\n").replace("\r", "\n").split("\n", -1);

        var anchors = new ArrayList<ASTNode>();
        collectAnchors(moduleNode, anchors);
        anchors.sort(Comparator.comparingLong(node -> start(node)));

        for( var range : ranges )
            applyRange(range, anchors, sourceLines, commentsMeta.getComments(), comments);
    }

    /**
     * Turn the `fmt: skip` / `fmt: off` / `fmt: on` comments into line ranges
     * excluded from formatting: `off` opens a range at its own line, closed by
     * the next `on` (or the end of the file if there is none); `skip` is a
     * one-line range.
     */
    private static List<int[]> computeRanges(List<Comment> allComments) {
        var ranges = new ArrayList<int[]>();
        int rangeStart = -1;
        for( var comment : allComments ) {
            var directive = fmtDirective(comment);
            if( directive == null )
                continue;
            if( rangeStart >= 0 ) {
                if( "on".equals(directive) ) {
                    ranges.add(new int[] { rangeStart, comment.line, 0 });
                    rangeStart = -1;
                }
                // a further "off" or "skip" while already excluded is a no-op
            }
            else if( "off".equals(directive) ) {
                rangeStart = comment.line;
            }
            else if( "skip".equals(directive) ) {
                ranges.add(new int[] { comment.line, comment.line, 1 });
            }
        }
        if( rangeStart >= 0 )
            ranges.add(new int[] { rangeStart, Integer.MAX_VALUE, 0 });
        return ranges;
    }

    private static String fmtDirective(Comment comment) {
        if( comment.line != comment.lastLine )
            return null;
        var matcher = FMT_DIRECTIVE.matcher(comment.text);
        return matcher.matches() ? matcher.group(1) : null;
    }

    /**
     * Resolve a single excluded range against the anchors (sorted in source
     * order): find the anchors it intersects, extend the range to cover them
     * completely, and assign the resulting verbatim text to the first of them
     * -- the rest are marked suppressed, since they emit nothing themselves.
     *
     * Every comment whose own line falls within the final range is marked as
     * emitted: it is printed as part of the raw verbatim text rather than via
     * the normal comment-printing path, so without this the `-format` safety
     * check ({@code getMissingComments()}) would see it as lost.
     */
    private void applyRange(int[] range, List<ASTNode> anchors, String[] sourceLines, List<Comment> allComments, CommentAttacher comments) {
        var startLine = range[0];
        var endLine = range[1];
        var skip = range[2] == 1;

        List<ASTNode> intersecting;
        if( skip ) {
            // `fmt: skip` applies to the outermost statement or declaration that
            // ENDS on the directive's line. A nested statement that merely shares
            // that last line (e.g. a closure body in `x = foo(bar { baz() }) // fmt: skip`)
            // must not steal the directive from the multi-line statement it is part of.
            ASTNode target = null;
            for( var anchor : anchors ) {
                if( anchor.getLineNumber() > endLine )
                    break;
                if( anchor.getLastLineNumber() == endLine && (target == null || start(anchor) < start(target)) )
                    target = anchor;
            }
            if( target == null )
                return;
            intersecting = List.of(target);
        }
        else {
            var overlapping = new ArrayList<ASTNode>();
            for( var anchor : anchors ) {
                // anchors are sorted by start, so once one begins past the range
                // no later one can intersect it either
                if( anchor.getLineNumber() > endLine )
                    break;
                if( anchor.getLastLineNumber() >= startLine )
                    overlapping.add(anchor);
            }

            // an anchor that merely encloses the range (e.g. the workflow around a
            // `fmt: off` region) is not part of it -- keep enclosing anchors only
            // when the range covers them entirely
            var maximal = new ArrayList<ASTNode>();
            for( var anchor : overlapping ) {
                var coveredByRange = anchor.getLineNumber() >= range[0] && anchor.getLastLineNumber() <= range[1];
                if( !coveredByRange && containsAnother(anchor, overlapping) )
                    continue;
                maximal.add(anchor);
            }
            intersecting = maximal;
        }
        if( intersecting.isEmpty() )
            return;

        for( var anchor : intersecting ) {
            startLine = Math.min(startLine, anchor.getLineNumber());
            endLine = Math.max(endLine, anchor.getLastLineNumber());
        }
        endLine = Math.min(endLine, sourceLines.length);
        // drop trailing empty lines (e.g. past the end of the file)
        while( endLine > startLine && sourceLines[endLine - 1].isEmpty() )
            endLine--;

        for( var comment : allComments ) {
            if( comment.line >= startLine && comment.line <= endLine ) {
                comments.markEmitted(comment);
                verbatimComments.add(comment);
            }
        }

        var text = String.join("\n", Arrays.copyOfRange(sourceLines, startLine - 1, endLine));
        var first = intersecting.get(0);
        verbatimSource.put(first, text);
        verbatimStart.put(first, startLine);
        for( int i = 1; i < intersecting.size(); i++ )
            suppressed.add(intersecting.get(i));
    }

    private static boolean containsAnother(ASTNode anchor, List<ASTNode> anchors) {
        for( var other : anchors ) {
            if( other == anchor )
                continue;
            if( start(anchor) <= start(other) && end(other) <= end(anchor) )
                return true;
        }
        return false;
    }

    /// ANCHOR COLLECTION
    //
    // An anchor is an AST node that the formatter emits through one of the
    // hooked emitters (see Formatter.appendVerbatim and its call sites): a
    // top-level script/config declaration, or a statement within a block body
    // (including a process/agent/config-apply directive line). This mirrors
    // CommentAttacher's traversal of containers and children, but only the
    // nodes that can themselves be printed verbatim need to be collected --
    // individual typed parameters, record fields and other finer-grained
    // constructs are out of scope for this PR (deferred, as noted in the spec).

    private final CodeVisitorSupport exprWalker = new CodeVisitorSupport() {
        @Override
        public void visitClosureExpression(ClosureExpression node) {
            collectStatements(node.getCode());
        }
    };

    private List<ASTNode> anchors;

    private void collectAnchors(ModuleNode moduleNode, List<ASTNode> anchors) {
        this.anchors = anchors;
        if( moduleNode instanceof ScriptNode sn )
            collectScript(sn);
        else if( moduleNode instanceof ConfigNode cn )
            collectConfigStatements(cn.getConfigStatements());
    }

    private void addAnchor(ASTNode node) {
        if( hasPosition(node) )
            anchors.add(node);
    }

    private void collectScript(ScriptNode node) {
        for( var decl : node.getDeclarations() ) {
            // a code snippet is an implicit entry workflow with no position
            if( decl instanceof WorkflowNode wn && wn.isCodeSnippet() ) {
                collectStatements(wn.main);
                continue;
            }
            if( !hasPosition(decl) )
                continue;
            anchors.add(decl);
            collectScriptDeclaration(decl);
        }
    }

    private void collectScriptDeclaration(ASTNode decl) {
        if( decl instanceof AgentNode an ) {
            collectDirectives(an.directives);
            collectStatements(an.prompt);
            // an.outputs uses typed-output formatting (not hooked) -- deferred
        }
        else if( decl instanceof FunctionNode fn ) {
            collectStatements(fn.getCode());
        }
        else if( decl instanceof OutputBlockNode obn ) {
            for( var on : obn.declarations )
                collectOutputBody(on.body);
        }
        else if( decl instanceof ProcessNodeV1 pn ) {
            collectDirectives(pn.directives);
            collectDirectives(pn.inputs);
            collectDirectives(pn.outputs);
            collectStatements(pn.exec);
            collectStatements(pn.stub);
        }
        else if( decl instanceof ProcessNodeV2 pn ) {
            collectDirectives(pn.directives);
            collectDirectives(pn.stagers);
            collectStatements(pn.exec);
            collectStatements(pn.stub);
            // pn.inputs/outputs/topics use typed formatting (not hooked) -- deferred
        }
        else if( decl instanceof WorkflowNode wn ) {
            collectStatements(wn.main);
            collectStatements(wn.onComplete);
            collectStatements(wn.onError);
            // wn.emits/publishers use typed-output formatting (not hooked) -- deferred
        }
        // ClassNode (enum), FeatureFlagNode, IncludeNode, ParamBlockNode,
        // ParamNodeV1 and RecordNode are leaf declarations with no nested
        // verbatim-hookable statements
    }

    /**
     * A block of directive-style statements (a bare method call per line,
     * e.g. a process directive), emitted through {@code Formatter.visitDirective}.
     */
    private void collectDirectives(Statement statement) {
        for( var stmt : asBlockStatements(statement) ) {
            if( asMethodCallX(stmt) != null )
                addAnchor(stmt);
        }
    }

    /**
     * The body of an `output` block declaration: a block of directives, except
     * for an `index { ... }` statement, whose own line is not individually
     * hookable but whose nested directives are.
     */
    private void collectOutputBody(Statement body) {
        for( var stmt : asBlockStatements(body) ) {
            var call = asMethodCallX(stmt);
            if( call == null )
                continue;
            if( "index".equals(call.getMethodAsString()) ) {
                var args = asMethodCallArguments(call);
                if( args.size() == 1 && args.get(0) instanceof ClosureExpression ce ) {
                    collectDirectives(ce.getCode());
                    continue;
                }
            }
            addAnchor(stmt);
        }
    }

    private void collectStatements(Statement node) {
        if( node == null )
            return;
        for( var stmt : asBlockStatements(node) )
            collectStatement(stmt);
    }

    private void collectStatement(Statement node) {
        if( node instanceof IfStatement is ) {
            addAnchor(is);
            collectStatements(is.getIfBlock());
            var elseBlock = is.getElseBlock();
            if( elseBlock instanceof IfStatement chained )
                collectElseIfChain(chained);
            else if( !(elseBlock instanceof EmptyStatement) )
                collectStatements(elseBlock);
            return;
        }
        if( node instanceof TryCatchStatement tcs ) {
            addAnchor(tcs);
            collectStatements(tcs.getTryStatement());
            for( var cs : tcs.getCatchStatements() )
                collectStatements(cs.getCode());
            return;
        }
        if( node instanceof ExpressionStatement || node instanceof ReturnStatement
                || node instanceof AssertStatement || node instanceof ThrowStatement ) {
            addAnchor(node);
        }
        // descend into the expression tree to reach any nested closure body
        node.visit(exprWalker);
    }

    /**
     * An `else if` link of a chain: unlike the root `if`, it is never
     * dispatched through {@code Formatter.visitIfElse} on its own (the parent
     * prints it inline), so it is not itself an anchor -- only the statements
     * in its branches are.
     */
    private void collectElseIfChain(IfStatement node) {
        collectStatements(node.getIfBlock());
        var elseBlock = node.getElseBlock();
        if( elseBlock instanceof IfStatement chained )
            collectElseIfChain(chained);
        else if( !(elseBlock instanceof EmptyStatement) )
            collectStatements(elseBlock);
    }

    /// CONFIG STATEMENTS

    private void collectConfigStatements(List<ConfigStatement> statements) {
        for( var stmt : statements ) {
            if( !hasPosition(stmt) )
                continue;
            anchors.add(stmt);
            if( stmt instanceof ConfigApplyBlockNode cabn ) {
                for( var apply : cabn.statements )
                    addAnchor(apply);
            }
            else if( stmt instanceof ConfigAssignNode can ) {
                can.value.visit(exprWalker);
            }
            else if( stmt instanceof ConfigBlockNode cbn ) {
                collectConfigStatements(cbn.statements);
            }
            else if( stmt instanceof ConfigIncludeNode cin ) {
                cin.source.visit(exprWalker);
            }
        }
    }

}
