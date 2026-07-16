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

import java.util.List;

import nextflow.script.ast.ASTNodeMarker;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.CodeVisitorSupport;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.Variable;
import org.codehaus.groovy.ast.expr.BinaryExpression;
import org.codehaus.groovy.ast.expr.BitwiseNegationExpression;
import org.codehaus.groovy.ast.expr.BooleanExpression;
import org.codehaus.groovy.ast.expr.CastExpression;
import org.codehaus.groovy.ast.expr.ClassExpression;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.ConstantExpression;
import org.codehaus.groovy.ast.expr.ConstructorCallExpression;
import org.codehaus.groovy.ast.expr.DeclarationExpression;
import org.codehaus.groovy.ast.expr.ElvisOperatorExpression;
import org.codehaus.groovy.ast.expr.EmptyExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.GStringExpression;
import org.codehaus.groovy.ast.expr.ListExpression;
import org.codehaus.groovy.ast.expr.MapEntryExpression;
import org.codehaus.groovy.ast.expr.MapExpression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.NamedArgumentListExpression;
import org.codehaus.groovy.ast.expr.NotExpression;
import org.codehaus.groovy.ast.expr.PropertyExpression;
import org.codehaus.groovy.ast.expr.RangeExpression;
import org.codehaus.groovy.ast.expr.TernaryExpression;
import org.codehaus.groovy.ast.expr.TupleExpression;
import org.codehaus.groovy.ast.expr.UnaryMinusExpression;
import org.codehaus.groovy.ast.expr.UnaryPlusExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.AssertStatement;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.CatchStatement;
import org.codehaus.groovy.ast.stmt.EmptyStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.IfStatement;
import org.codehaus.groovy.ast.stmt.ReturnStatement;
import org.codehaus.groovy.ast.stmt.ThrowStatement;
import org.codehaus.groovy.ast.stmt.TryCatchStatement;
import org.codehaus.groovy.runtime.DefaultGroovyMethods;
import org.codehaus.groovy.syntax.Types;

import static nextflow.script.ast.ASTUtils.*;

/**
 * Transform an AST into formatted source code.
 *
 * Formatting never removes comments (expects comment metadata prepared
 * by {@link CommentReattacher}):
 *
 * <ul>
 *   <li>trailing comments are emitted on the same line as their statement;</li>
 *   <li>dangling comments (after the last statement of a block, or in an
 *       empty block) are emitted before the closing brace;</li>
 *   <li>closures with comments are always formatted as multi-line blocks
 *       so that their comments can be emitted;</li>
 *   <li>leading comments are normalized against Windows (CRLF) line
 *       endings.</li>
 * </ul>
 *
 * @see CommentReattacher
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class Formatter extends CodeVisitorSupport {

    private FormattingOptions options;

    private StringBuilder builder = new StringBuilder();

    private int indentCount = 0;

    public Formatter(FormattingOptions options) {
        this.options = options;
    }

    public void append(char c) {
        builder.append(c);
    }

    public void append(String str) {
        builder.append(str);
    }

    public void appendIndent() {
        var str = options.insertSpaces()
            ? " ".repeat(options.tabSize() * indentCount)
            : "\t".repeat(indentCount);
        builder.append(str);
    }

    public void appendNewLine() {
        builder.append('\n');
    }

    public void appendLeadingComments(ASTNode node) {
        var comments = (List<String>) node.getNodeMetaData(ASTNodeMarker.LEADING_COMMENTS);
        if( comments == null || comments.isEmpty() )
            return;

        for( var line : DefaultGroovyMethods.asReversed(comments) ) {
            if( isNewLine(line) ) {
                appendNewLine();
            }
            else {
                appendIndent();
                append(stripCarriageReturns(line.stripLeading()));
            }
        }
    }

    public boolean hasLeadingComments(ASTNode node) {
        var comments = (List<String>) node.getNodeMetaData(ASTNodeMarker.LEADING_COMMENTS);
        return comments != null && !comments.isEmpty();
    }

    public boolean hasTrailingComment(ASTNode node) {
        var comment = (String) node.getNodeMetaData(ASTNodeMarker.TRAILING_COMMENT);
        return comment != null;
    }

    public void appendTrailingComment(ASTNode node) {
        var comment = (String) node.getNodeMetaData(ASTNodeMarker.TRAILING_COMMENT);
        if( comment != null ) {
            append(' ');
            append(stripCarriageReturns(comment));
        }
    }

    /**
     * Emit the comments that follow the last statement of a block (or fill
     * an otherwise-empty block), before the closing brace.
     */
    public void appendDanglingComments(ASTNode node) {
        appendCommentList((List<String>) node.getNodeMetaData(CommentReattacher.DANGLING_COMMENTS));
    }

    /**
     * Emit the comments that follow a node at the end of a section that has
     * no block emission point (e.g. the last workflow take).
     */
    public void appendDanglingAfter(ASTNode node) {
        appendCommentList((List<String>) node.getNodeMetaData(CommentReattacher.DANGLING_AFTER));
    }

    private void appendCommentList(List<String> comments) {
        if( comments == null || comments.isEmpty() )
            return;
        var pendingNewLine = false;
        for( var line : comments ) {
            if( isNewLine(line) ) {
                appendNewLine();
                pendingNewLine = false;
            }
            else {
                if( pendingNewLine )
                    appendNewLine();
                appendIndent();
                append(stripCarriageReturns(line.stripLeading()));
                pendingNewLine = true;
            }
        }
        if( pendingNewLine )
            appendNewLine();
    }

    public boolean hasDanglingComments(ASTNode node) {
        var comments = (List<String>) node.getNodeMetaData(CommentReattacher.DANGLING_COMMENTS);
        return comments != null && !comments.isEmpty();
    }

    /**
     * Emit a node excluded from formatting with `fmt: skip` or `fmt: off`
     * verbatim from the source text (or nothing, if the node is part of a
     * verbatim region emitted by an earlier node). Returns false if the
     * node is not part of a verbatim region.
     */
    public boolean appendVerbatim(ASTNode node) {
        if( node.getNodeMetaData(CommentReattacher.VERBATIM_SUPPRESSED) != null )
            return true;
        var text = (String) node.getNodeMetaData(CommentReattacher.VERBATIM_SOURCE);
        if( text == null )
            return false;
        appendLeadingComments(node);
        append(text);
        appendNewLine();
        return true;
    }

    private static boolean isNewLine(String line) {
        return "\n".equals(line) || "\r\n".equals(line);
    }

    private static String stripCarriageReturns(String line) {
        return line.replace("\r", "");
    }

    public void incIndent() {
        indentCount++;
    }

    public void decIndent() {
        indentCount--;
    }

    public String toString() {
        return builder.toString();
    }

    // line-length wrapping

    /**
     * How aggressively expressions are being wrapped onto multiple lines,
     * on top of preserving the wrapping in the source:
     *
     * NONE -- only expressions that were wrapped in the source are wrapped;
     * TOP  -- the outermost wrappable construct of the statement is wrapped;
     * ALL  -- every wrappable construct of the statement is wrapped.
     */
    private enum WrapMode {
        NONE,
        TOP,
        ALL
    }

    private WrapMode wrapMode = WrapMode.NONE;

    private int exprDepth = 0;

    /**
     * Emit a line of output, and if it exceeds the maximum line length,
     * roll it back and re-emit it with expressions wrapped onto multiple
     * lines -- first only the outermost construct, then, if the result
     * still has an over-long line (e.g. deeply nested arguments), every
     * construct.
     */
    public void emitWrappable(Runnable body) {
        if( options.maxLineLength() <= 0 || wrapMode != WrapMode.NONE ) {
            body.run();
            return;
        }
        var mark = builder.length();
        body.run();
        if( !exceedsMaxLineLength(mark) )
            return;
        builder.setLength(mark);
        wrapMode = WrapMode.TOP;
        body.run();
        if( exceedsMaxLineLength(mark) ) {
            builder.setLength(mark);
            wrapMode = WrapMode.ALL;
            body.run();
        }
        wrapMode = WrapMode.NONE;
    }

    /**
     * Determine whether any line emitted after the given position exceeds
     * the maximum line length.
     */
    private boolean exceedsMaxLineLength(int from) {
        var max = options.maxLineLength();
        var lineStart = builder.lastIndexOf("\n", from - 1) + 1;
        var col = 0;
        for( int i = lineStart; i < builder.length(); i++ ) {
            var c = builder.charAt(i);
            if( c == '\n' )
                col = 0;
            else if( c == '\t' )
                col += options.tabSize();
            else
                col++;
            if( col > max )
                return true;
        }
        return false;
    }

    /**
     * Determine whether the current construct should be wrapped because of
     * the line-length limit.
     */
    private boolean forceWrap() {
        if( wrapMode == WrapMode.ALL )
            return true;
        if( wrapMode == WrapMode.TOP )
            return exprDepth <= 2;
        return false;
    }

    // statements

    @Override
    public void visitBlockStatement(BlockStatement node) {
        super.visitBlockStatement(node);
        appendDanglingComments(node);
    }

    @Override
    public void visitIfElse(IfStatement node) {
        if( appendVerbatim(node) )
            return;
        appendLeadingComments(node);
        appendIndent();
        visitIfElse(node, true);
        appendTrailingComment(node);
        appendNewLine();
    }

    /**
     * Format an if/else statement in K&R style ("} else {"), without the
     * line terminator (the caller appends the trailing comment and newline
     * after the final closing brace).
     */
    protected void visitIfElse(IfStatement node, boolean isFirst) {
        append("if ");
        visit(node.getBooleanExpression());
        append(" {\n");
        incIndent();
        visit(node.getIfBlock());
        decIndent();
        appendIndent();
        append('}');
        if( node.getElseBlock() instanceof IfStatement is ) {
            append(" else ");
            visitIfElse(is, false);
        }
        else if( !(node.getElseBlock() instanceof EmptyStatement) ) {
            append(" else {\n");
            incIndent();
            visit(node.getElseBlock());
            decIndent();
            appendIndent();
            append('}');
        }
    }

    private Expression currentRootExpr;

    @Override
    public void visitExpressionStatement(ExpressionStatement node) {
        if( appendVerbatim(node) )
            return;
        appendLeadingComments(node);
        emitWrappable(() -> {
            var cre = currentRootExpr;
            currentRootExpr = node.getExpression();
            appendIndent();
            visitStatementLabels(node);
            visit(node.getExpression());
            appendTrailingComment(node);
            appendNewLine();
            currentRootExpr = cre;
        });
    }

    private void visitStatementLabels(ExpressionStatement node) {
        if( node.getStatementLabels() == null )
            return;
        for( var label : DefaultGroovyMethods.asReversed(node.getStatementLabels()) ) {
            append(label);
            append(": ");
        }
    }

    @Override
    public void visitReturnStatement(ReturnStatement node) {
        if( appendVerbatim(node) )
            return;
        appendLeadingComments(node);
        emitWrappable(() -> {
            var cre = currentRootExpr;
            currentRootExpr = node.getExpression();
            appendIndent();
            append("return ");
            visit(node.getExpression());
            appendTrailingComment(node);
            appendNewLine();
            currentRootExpr = cre;
        });
    }

    @Override
    public void visitAssertStatement(AssertStatement node) {
        if( appendVerbatim(node) )
            return;
        appendLeadingComments(node);
        emitWrappable(() -> {
            appendIndent();
            append("assert ");
            visit(node.getBooleanExpression());
            if( !(node.getMessageExpression() instanceof ConstantExpression ce && ce.isNullExpression()) ) {
                append(" : ");
                visit(node.getMessageExpression());
            }
            appendTrailingComment(node);
            appendNewLine();
        });
    }

    @Override
    public void visitTryCatchFinally(TryCatchStatement node) {
        if( appendVerbatim(node) )
            return;
        appendLeadingComments(node);
        appendIndent();
        append("try {\n");
        incIndent();
        visit(node.getTryStatement());
        decIndent();
        appendIndent();
        append('}');
        for( var catchStatement : node.getCatchStatements() ) {
            visitCatchStatement(catchStatement);
        }
        appendTrailingComment(node);
        appendNewLine();
    }

    @Override
    public void visitThrowStatement(ThrowStatement node) {
        if( appendVerbatim(node) )
            return;
        appendLeadingComments(node);
        emitWrappable(() -> {
            appendIndent();
            append("throw ");
            visit(node.getExpression());
            appendTrailingComment(node);
            appendNewLine();
        });
    }

    @Override
    public void visitCatchStatement(CatchStatement node) {
        // format in K&R style ("} catch (...) {"), unless the catch clause
        // has leading comments that must be emitted on their own lines
        if( node.getNodeMetaData(ASTNodeMarker.LEADING_COMMENTS) != null ) {
            appendNewLine();
            appendLeadingComments(node);
            appendIndent();
        }
        else {
            append(' ');
        }
        append("catch (");

        var variable = node.getVariable();
        append(variable.getName());
        if( hasType(variable) ) {
            append(": ");
            append(variable.getType().getNameWithoutPackage());
        }

        append(") {\n");
        incIndent();
        visit(node.getCode());
        decIndent();
        appendIndent();
        append('}');
    }

    // expressions

    private boolean inWrappedMethodChain;

    @Override
    public void visitMethodCallExpression(MethodCallExpression node) {
        var beginWrappedMethodChain = shouldWrapMethodChain(node);
        if( beginWrappedMethodChain )
            inWrappedMethodChain = true;

        if( !node.isImplicitThis() ) {
            var receiver = node.getObjectExpression();
            visit(receiver);
            if( inWrappedMethodChain ) {
                incIndent();
                if( !nextflow.script.dsl.Types.isNamespace(receiver.getType()) || hasLeadingComments(node) ) {
                    appendNewLine();
                    appendLeadingComments(node);
                    appendIndent();
                }
            }
            if( node.isSpreadSafe() )
                append('*');
            else if( node.isSafe() )
                append('?');
            append('.');
        }

        visit(node.getMethod());

        var iwmc = inWrappedMethodChain;
        inWrappedMethodChain = false;
        var args = asMethodCallArguments(node);
        var lastClosureArg = args.size() > 0 && args.get(args.size() - 1) instanceof ClosureExpression;
        var parenArgs = lastClosureArg
            ? DefaultGroovyMethods.init(args)
            : args;
        if( parenArgs.size() > 0 || !lastClosureArg ) {
            var wrap = shouldWrapMethodCall(node);
            append('(');
            if( wrap )
                incIndent();
            visitArguments(parenArgs, wrap);
            if( wrap ) {
                appendNewLine();
                decIndent();
                appendIndent();
            }
            append(')');
        }
        if( lastClosureArg ) {
            append(' ');
            visit(args.get(args.size() - 1));
        }
        inWrappedMethodChain = iwmc;

        if( !node.isImplicitThis() && inWrappedMethodChain )
            decIndent();
        if( beginWrappedMethodChain )
            inWrappedMethodChain = false;
    }

    @Override
    public void visitConstructorCallExpression(ConstructorCallExpression node) {
        append("new ");
        visitTypeAnnotation(node.getType());
        append('(');
        visitArguments(asMethodCallArguments(node), false);
        append(')');
    }

    public void visitDirective(MethodCallExpression call) {
        visitDirective(call, null);
    }

    /**
     * Format a directive, emitting the trailing comment of the given source
     * node (typically the enclosing statement) at the end of the line.
     */
    public void visitDirective(MethodCallExpression call, ASTNode trailingSource) {
        appendIndent();
        append(call.getMethodAsString());
        var arguments = asMethodCallArguments(call);
        if( !arguments.isEmpty() ) {
            append(' ');
            visitArguments(arguments, false);
        }
        if( trailingSource != null )
            appendTrailingComment(trailingSource);
        appendNewLine();
    }

    public void visitArguments(List<Expression> args, boolean wrap) {
        var hasNamedArgs = args.size() > 0 && args.get(0) instanceof NamedArgumentListExpression;
        var positionalArgs = hasNamedArgs
            ? DefaultGroovyMethods.tail(args)
            : args;
        visitPositionalArgs(positionalArgs, wrap);
        if( hasNamedArgs ) {
            if( positionalArgs.size() > 0 )
                append(wrap ? "," : ", ");
            var mapX = (MapExpression)args.get(0);
            visitNamedArgs(mapX.getMapEntryExpressions(), wrap);
        }
    }

    private boolean inVariableDeclaration;

    @Override
    public void visitBinaryExpression(BinaryExpression node) {
        if( node instanceof DeclarationExpression ) {
            append("def ");
            inVariableDeclaration = true;
            visit(node.getLeftExpression());
            inVariableDeclaration = false;
            var source = node.getRightExpression();
            if( !(source instanceof EmptyExpression) ) {
                append(" = ");
                var cre = currentRootExpr;
                currentRootExpr = source;
                visit(source);
                currentRootExpr = cre;
            }
            return;
        }

        if( node.getOperation().isA(Types.LEFT_SQUARE_BRACKET) ) {
            visit(node.getLeftExpression());
            append('[');
            visit(node.getRightExpression());
            append(']');
            return;
        }

        visit(node.getLeftExpression());

        append(' ');
        append(node.getOperation().getText());
        append(' ');

        if( node.getOperation().isA(Types.ASSIGNMENT_OPERATOR) ) {
            var source = node.getRightExpression();
            var cre = currentRootExpr;
            currentRootExpr = source;
            visit(source);
            currentRootExpr = cre;
        }
        else {
            visit(node.getRightExpression());
        }
    }

    @Override
    public void visitTernaryExpression(TernaryExpression node) {
        if( shouldWrapExpression(node) ) {
            visit(node.getBooleanExpression());
            incIndent();
            appendNewLine();
            appendIndent();
            append("? ");
            visit(node.getTrueExpression());
            appendNewLine();
            appendIndent();
            append(": ");
            visit(node.getFalseExpression());
            decIndent();
        }
        else {
            visit(node.getBooleanExpression());
            append(" ? ");
            visit(node.getTrueExpression());
            append(" : ");
            visit(node.getFalseExpression());
        }
    }

    @Override
    public void visitShortTernaryExpression(ElvisOperatorExpression node) {
        visit(node.getTrueExpression());
        append(" ?: ");
        visit(node.getFalseExpression());
    }

    @Override
    public void visitBooleanExpression(BooleanExpression node) {
        visit(node.getExpression());
    }

    @Override
    public void visitNotExpression(NotExpression node) {
        append('!');
        visit(node.getExpression());
    }

    @Override
    public void visitClosureExpression(ClosureExpression node) {
        append('{');
        if( node.getParameters() != null && node.getParameters().length > 0 ) {
            append(' ');
            visitParameters(node.getParameters());
            append(" ->");
        }
        var code = (BlockStatement) node.getCode();
        if( code.getStatements().size() == 0 && !hasDanglingComments(code) ) {
            append(" }");
        }
        else if( code.getStatements().size() == 1
                && code.getStatements().get(0) instanceof ExpressionStatement es
                && !shouldWrapExpression(node)
                && !hasComments(code, es) ) {
            append(' ');
            visitStatementLabels(es);
            visit(es.getExpression());
            append(" }");
        }
        else {
            appendNewLine();
            incIndent();
            visit(code);
            decIndent();
            appendIndent();
            append('}');
        }
    }

    /**
     * Determine whether a single-statement closure carries comments and
     * must therefore be formatted as a multi-line block so that the
     * comments can be emitted.
     */
    private boolean hasComments(BlockStatement code, ExpressionStatement es) {
        return hasDanglingComments(code)
            || es.getNodeMetaData(ASTNodeMarker.LEADING_COMMENTS) != null
            || hasTrailingComment(es);
    }

    public void visitParameters(Parameter[] parameters) {
        for( int i = 0; i < parameters.length; i++ ) {
            var param = parameters[i];
            append(param.getName());
            if( hasType(param) ) {
                append(": ");
                visitTypeAnnotation(param.getType());
            }
            if( param.hasInitialExpression() ) {
                append(" = ");
                visit(param.getInitialExpression());
            }
            if( i + 1 < parameters.length )
                append(", ");
        }
    }

    @Override
    public void visitTupleExpression(TupleExpression node) {
        var wrap = hasTrailingComma(node) || (forceWrap() && node.getExpressions().size() > 1) || hasElementComments(node.getExpressions());
        append('(');
        if( wrap )
            incIndent();
        visitPositionalArgs(node.getExpressions(), wrap);
        if( wrap ) {
            appendNewLine();
            decIndent();
            appendIndent();
        }
        append(')');
    }

    @Override
    public void visitListExpression(ListExpression node) {
        var wrap = hasTrailingComma(node) || (forceWrap() && node.getExpressions().size() > 1) || hasElementComments(node.getExpressions());
        append('[');
        if( wrap )
            incIndent();
        visitPositionalArgs(node.getExpressions(), wrap);
        if( wrap ) {
            appendNewLine();
            decIndent();
            appendIndent();
        }
        append(']');
    }

    protected void visitPositionalArgs(List<Expression> args, boolean wrap) {
        var comma = wrap ? "," : ", ";
        var trailingComma = wrap && args.size() > 1;
        for( int i = 0; i < args.size(); i++ ) {
            if( wrap ) {
                appendNewLine();
                appendLeadingComments(args.get(i));
                appendIndent();
            }
            visit(args.get(i));
            if( trailingComma || i + 1 < args.size() )
                append(comma);
        }
    }

    @Override
    public void visitMapExpression(MapExpression node) {
        if( node.getMapEntryExpressions().isEmpty() ) {
            append("[:]");
            return;
        }
        var wrap = hasTrailingComma(node) || (forceWrap() && node.getMapEntryExpressions().size() > 1) || hasElementComments(node.getMapEntryExpressions());
        append('[');
        if( wrap )
            incIndent();
        visitNamedArgs(node.getMapEntryExpressions(), wrap);
        if( wrap ) {
            appendNewLine();
            decIndent();
            appendIndent();
        }
        append(']');
    }

    protected void visitNamedArgs(List<MapEntryExpression> args, boolean wrap) {
        var comma = wrap ? "," : ", ";
        var trailingComma = wrap && args.size() > 1;
        for( int i = 0; i < args.size(); i++ ) {
            if( wrap ) {
                appendNewLine();
                appendLeadingComments(args.get(i));
                appendIndent();
            }
            visit(args.get(i));
            if( trailingComma || i + 1 < args.size() )
                append(comma);
        }
    }

    @Override
    public void visitMapEntryExpression(MapEntryExpression node) {
        visit(node.getKeyExpression());
        append(": ");
        visit(node.getValueExpression());
    }

    @Override
    public void visitRangeExpression(RangeExpression node) {
        visit(node.getFrom());
        if( node.isExclusiveLeft() )
            append('<');
        append("..");
        if( node.isExclusiveRight() )
            append('<');
        visit(node.getTo());
    }

    @Override
    public void visitUnaryMinusExpression(UnaryMinusExpression node) {
        append('-');
        visit(node.getExpression());
    }

    @Override
    public void visitUnaryPlusExpression(UnaryPlusExpression node) {
        append('+');
        visit(node.getExpression());
    }

    @Override
    public void visitBitwiseNegationExpression(BitwiseNegationExpression node) {
        append('~');
        visit(node.getExpression());
    }

    @Override
    public void visitCastExpression(CastExpression node) {
        visit(node.getExpression());
        append(" as ");
        visitTypeAnnotation(node.getType());
    }

    @Override
    public void visitConstantExpression(ConstantExpression node) {
        var text = (String) node.getNodeMetaData(ASTNodeMarker.VERBATIM_TEXT);
        if( text != null )
            append(reindentMultiLineString(text, indentDelta(node)));
        else
            append(node.getText());
    }

    /**
     * The current column (0-based) at the end of the emitted text.
     */
    private int currentColumn() {
        return builder.length() - builder.lastIndexOf("\n") - 1;
    }

    /**
     * The difference between the column where a multi-line string is being
     * emitted and the column where it appeared in the source, so that its
     * interior lines can be shifted accordingly.
     */
    private int indentDelta(org.codehaus.groovy.ast.expr.Expression node) {
        if( !options.insertSpaces() )
            return 0;
        if( node.getLineNumber() < 0 || node.getColumnNumber() < 0 )
            return 0;
        return currentColumn() - (node.getColumnNumber() - 1);
    }

    /**
     * Shift the interior lines of a multi-line string so that its
     * indentation follows the indentation of the emitted code. Blank lines
     * are left untouched, and lines are never shifted past their first
     * non-whitespace character.
     */
    private static String reindentMultiLineString(String text, int delta) {
        if( delta == 0 || text.indexOf('\n') < 0 )
            return text;
        var lines = text.split("\n", -1);
        var result = new StringBuilder(lines[0]);
        for( int i = 1; i < lines.length; i++ ) {
            result.append('\n');
            var line = lines[i];
            if( line.isEmpty() ) {
                // empty lines are left untouched; whitespace-only lines are
                // shifted since they carry the indentation of whatever
                // follows (an interpolation or the closing quote)
            }
            else if( delta > 0 ) {
                result.append(" ".repeat(delta));
                result.append(line);
            }
            else {
                int strip = 0;
                while( strip < -delta && strip < line.length() && line.charAt(strip) == ' ' )
                    strip++;
                result.append(line.substring(strip));
            }
        }
        return result.toString();
    }

    @Override
    public void visitClassExpression(ClassExpression node) {
        visitTypeAnnotation(node.getType());
    }

    public void visitTypeAnnotation(ClassNode type) {
        if( isLegacyType(type) )
            append(type.getNodeMetaData(ASTNodeMarker.LEGACY_TYPE));
        else
            append(nextflow.script.dsl.Types.getName(type));
    }

    @Override
    public void visitVariableExpression(VariableExpression node) {
        append(node.getText());
        if( inVariableDeclaration && hasType(node) ) {
            append(": ");
            visitTypeAnnotation(node.getType());
        }
    }

    @Override
    public void visitPropertyExpression(PropertyExpression node) {
        visit(node.getObjectExpression());
        if( node.isSpreadSafe() )
            append('*');
        else if( node.isSafe() )
            append('?');
        append('.');
        visit(node.getProperty());
    }

    @Override
    public void visitGStringExpression(GStringExpression node) {
        // see also: GStringUtil.writeToImpl()
        var quoteChar = (String) node.getNodeMetaData(ASTNodeMarker.QUOTE_CHAR, k -> DQ_STR);
        var delta = indentDelta(node);
        append(quoteChar);
        var ss = node.getStrings();
        var vs = node.getValues();
        for( int i = 0; i < ss.size(); i++ ) {
            var string = ss.get(i);
            var text = (String) string.getNodeMetaData(ASTNodeMarker.VERBATIM_TEXT);
            if( text != null )
                append(reindentMultiLineString(text, delta));
            if( i < vs.size() ) {
                append("${");
                visit(vs.get(i));
                append('}');
            }
        }
        append(quoteChar);
    }

    @Override
    public void visit(Expression node) {
        var number = (Number) node.getNodeMetaData(ASTNodeMarker.INSIDE_PARENTHESES_LEVEL);
        if( number != null && number.intValue() > 0 )
            append('(');
        exprDepth++;
        super.visit(node);
        exprDepth--;
        if( number != null && number.intValue() > 0 )
            append(')');
    }

    // helpers

    private static boolean hasTrailingComma(Expression node) {
        return node.getNodeMetaData(ASTNodeMarker.TRAILING_COMMA) != null;
    }

    public static boolean hasType(ClassNode type) {
        return !ClassHelper.isDynamicTyped(type) || isLegacyType(type);
    }

    public static boolean hasType(Variable variable) {
        return !variable.isDynamicTyped() || isLegacyType(variable.getType());
    }

    public static boolean isLegacyType(ClassNode cn) {
        return cn.getNodeMetaData(ASTNodeMarker.LEGACY_TYPE) != null;
    }

    private boolean shouldWrapExpression(Expression node) {
        return node.getLineNumber() < node.getLastLineNumber() || forceWrap();
    }

    private boolean shouldWrapMethodCall(MethodCallExpression node) {
        if( hasTrailingComma(node.getArguments()) )
            return true;
        if( forceWrap() && asMethodCallArguments(node).size() > 0 )
            return true;
        for( var arg : asMethodCallArguments(node) ) {
            if( hasLeadingComments(arg) )
                return true;
            if( arg instanceof MapExpression named && hasElementComments(named.getMapEntryExpressions()) )
                return true;
        }
        var start = node.getMethod();
        var end = node.getArguments();
        return start.getLineNumber() < end.getLastLineNumber();
    }

    private boolean shouldWrapMethodChain(MethodCallExpression node) {
        if( currentRootExpr != node )
            return false;
        if( chainHasComments(node) )
            return true;
        if( !shouldWrapExpression(node) )
            return false;

        Expression root = node;
        int depth = 0;
        while( root instanceof MethodCallExpression mce && !mce.isImplicitThis() ) {
            root = mce.getObjectExpression();
            depth += 1;
        }

        if( wrapMode != WrapMode.NONE )
            return depth >= 2;

        return shouldWrapExpression(root)
            ? false
            : depth >= 2;
    }

    /**
     * Determine whether any link of a method chain carries comments (which
     * can only be emitted when the chain is wrapped).
     */
    private boolean chainHasComments(MethodCallExpression node) {
        Expression expr = node;
        while( expr instanceof MethodCallExpression mce && !mce.isImplicitThis() ) {
            if( hasLeadingComments(mce) )
                return true;
            expr = mce.getObjectExpression();
        }
        return false;
    }

    /**
     * Determine whether any element of a wrapped construct carries comments
     * (which can only be emitted when the construct is wrapped).
     */
    private boolean hasElementComments(List<? extends ASTNode> elements) {
        for( var element : elements ) {
            if( hasLeadingComments(element) )
                return true;
        }
        return false;
    }

    private static final String SLASH_STR = "/";
    private static final String TDQ_STR = "\"\"\"";
    private static final String TSQ_STR = "'''";
    private static final String SQ_STR = "'";
    private static final String DQ_STR = "\"";

}
