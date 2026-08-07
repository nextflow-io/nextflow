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
import java.util.stream.Stream;

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
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class Formatter extends CodeVisitorSupport {

    private FormattingOptions options;

    private StringBuilder builder = new StringBuilder();

    private int indentCount = 0;

    private CommentAttacher comments = CommentAttacher.of(null);

    private int stringIndentDelta = 0;

    public Formatter(FormattingOptions options) {
        this.options = options;
    }

    public void setComments(CommentAttacher comments) {
        this.comments = comments;
    }

    public CommentAttacher getComments() {
        return comments;
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

    /**
     * Append the comments that precede a node, along with any blank lines
     * that separate them from each other and from the node.
     *
     * @param node
     */
    public void appendLeadingComments(ASTNode node) {
        for( var comment : comments.leading(node) ) {
            appendBlankLines(comment.line);
            appendComment(comment);
        }
        if( node != null )
            appendBlankLines(node.getLineNumber());
    }

    /**
     * Append the comments that belong to a construct but not to any
     * particular child, e.g. a comment before a closing brace.
     *
     * @param node
     */
    public void appendDanglingComments(ASTNode node) {
        for( var comment : comments.dangling(node) ) {
            appendBlankLines(comment.line);
            appendComment(comment);
        }
    }

    public boolean hasTrailingComment(ASTNode node) {
        return !comments.trailing(node).isEmpty();
    }

    public void appendTrailingComment(ASTNode node) {
        for( var comment : comments.trailing(node) ) {
            this.comments.markEmitted(comment);
            append(' ');
            append(comment.text.replace('\n', ' '));
        }
    }

    private void appendBlankLines(int line) {
        if( builder.length() == 0 )
            return;
        var count = comments.blankLinesBefore(line);
        for( int i = 0; i < count; i++ )
            appendNewLine();
    }

    private void appendComment(Comment comment) {
        comments.markEmitted(comment);
        var lines = comment.text.split("\n", -1);
        for( int i = 0; i < lines.length; i++ ) {
            var line = lines[i].strip();
            appendIndent();
            // preserve the canonical alignment of a block comment
            if( i > 0 && line.startsWith("*") )
                append(' ');
            append(line);
            appendNewLine();
        }
    }

    /**
     * Begin emitting a statement, recording how far its indentation moved
     * from the source. A multi-line string inside the statement is shifted by
     * the same amount, so that the string body keeps its position relative to
     * the code around it -- which is the property that matters, since a
     * process script is dedented at runtime.
     *
     * Only a statement that begins its own source line has an indentation
     * to compare; for one that starts mid-line the column is a position within
     * a line of code, so its strings are left alone.
     *
     * @param node
     * @return the enclosing statement's shift, to be passed to {@link #endStatement}
     */
    public int beginStatement(ASTNode node) {
        var previous = stringIndentDelta;
        stringIndentDelta = options.insertSpaces() && comments.startsLine(node)
            ? options.tabSize() * indentCount - (node.getColumnNumber() - 1)
            : 0;
        return previous;
    }

    public void endStatement(int previous) {
        stringIndentDelta = previous;
    }

    /**
     * Shift the interior lines of a multi-line string by the indentation
     * shift of the statement that contains it. Blank lines are left alone,
     * and a line is never shifted past its first non-whitespace character.
     */
    private String reindentString(String text) {
        var delta = stringIndentDelta;
        if( delta == 0 || text.indexOf('\n') < 0 )
            return text;

        var lines = text.split("\n", -1);
        var builder = new StringBuilder(lines[0]);
        for( int i = 1; i < lines.length; i++ ) {
            builder.append('\n');
            var line = lines[i];
            if( line.isEmpty() ) {
                // an empty line has no indentation to shift
            }
            else if( delta > 0 ) {
                builder.append(" ".repeat(delta));
                builder.append(line);
            }
            else {
                int strip = 0;
                while( strip < -delta && strip < line.length() && line.charAt(strip) == ' ' )
                    strip++;
                builder.append(line, strip, line.length());
            }
        }
        return builder.toString();
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

    // statements

    @Override
    public void visitBlockStatement(BlockStatement node) {
        super.visitBlockStatement(node);
        appendDanglingComments(node);
    }

    @Override
    public void visitIfElse(IfStatement node) {
        visitIfElse(node, true);
    }

    protected void visitIfElse(IfStatement node, boolean preIndent) {
        appendLeadingComments(node);
        if( preIndent )
            appendIndent();
        append("if ");
        visit(node.getBooleanExpression());
        append(" {\n");
        incIndent();
        visit(node.getIfBlock());
        decIndent();
        appendIndent();
        append("}");
        var elseBlock = node.getElseBlock();
        if( elseBlock instanceof IfStatement is ) {
            appendTrailingComment(node);
            appendNewLine();
            appendIndent();
            append("else ");
            visitIfElse(is, false);
            return;
        }
        if( !(elseBlock instanceof EmptyStatement) ) {
            appendNewLine();
            appendIndent();
            append("else {\n");
            incIndent();
            visit(elseBlock);
            decIndent();
            appendIndent();
            append("}");
        }
        appendTrailingComment(node);
        appendNewLine();
    }

    private Expression currentRootExpr;

    @Override
    public void visitExpressionStatement(ExpressionStatement node) {
        var cre = currentRootExpr;
        currentRootExpr = node.getExpression();
        appendLeadingComments(node);
        appendIndent();
        var sid = beginStatement(node);
        visitStatementLabels(node);
        visit(node.getExpression());
        endStatement(sid);
        appendTrailingComment(node);
        appendNewLine();
        currentRootExpr = cre;
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
        var cre = currentRootExpr;
        currentRootExpr = node.getExpression();
        appendLeadingComments(node);
        appendIndent();
        var sid = beginStatement(node);
        append("return ");
        visit(node.getExpression());
        endStatement(sid);
        appendTrailingComment(node);
        appendNewLine();
        currentRootExpr = cre;
    }

    @Override
    public void visitAssertStatement(AssertStatement node) {
        appendLeadingComments(node);
        appendIndent();
        var sid = beginStatement(node);
        append("assert ");
        visit(node.getBooleanExpression());
        if( !(node.getMessageExpression() instanceof ConstantExpression ce && ce.isNullExpression()) ) {
            append(" : ");
            visit(node.getMessageExpression());
        }
        endStatement(sid);
        appendTrailingComment(node);
        appendNewLine();
    }

    @Override
    public void visitTryCatchFinally(TryCatchStatement node) {
        appendLeadingComments(node);
        appendIndent();
        append("try {\n");
        incIndent();
        visit(node.getTryStatement());
        decIndent();
        appendIndent();
        append("}");
        appendTrailingComment(node);
        appendNewLine();
        for( var catchStatement : node.getCatchStatements() ) {
            visit(catchStatement);
        }
    }

    @Override
    public void visitThrowStatement(ThrowStatement node) {
        appendLeadingComments(node);
        appendIndent();
        var sid = beginStatement(node);
        append("throw ");
        visit(node.getExpression());
        endStatement(sid);
        appendTrailingComment(node);
        appendNewLine();
    }

    @Override
    public void visitCatchStatement(CatchStatement node) {
        appendLeadingComments(node);
        appendIndent();
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
        append("}");
        appendTrailingComment(node);
        appendNewLine();
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
                var receiverComment = hasTrailingComment(receiver);
                appendTrailingComment(receiver);
                incIndent();
                var linkComments = comments.leading(node);
                // a namespace receiver stays on the same line, unless a
                // comment needs a line of its own
                if( !nextflow.script.dsl.Types.isNamespace(receiver.getType()) || !linkComments.isEmpty() || receiverComment ) {
                    appendNewLine();
                    for( var comment : linkComments )
                        appendComment(comment);
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
            var wrap = shouldWrapMethodCall(node) || hasComments(node, parenArgs);
            append('(');
            if( wrap )
                incIndent();
            visitArguments(parenArgs, wrap);
            if( wrap ) {
                appendNewLine();
                appendDanglingComments(node);
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
        visitDirective(call, call);
    }

    public void visitDirective(ASTNode owner, MethodCallExpression call) {
        appendIndent();
        var sid = beginStatement(owner);
        append(call.getMethodAsString());
        var arguments = asMethodCallArguments(call);
        if( !arguments.isEmpty() ) {
            var wrap = hasComments(call, arguments);
            if( wrap )
                incIndent();
            else
                append(' ');
            visitArguments(arguments, wrap);
            if( wrap ) {
                appendNewLine();
                appendDanglingComments(call);
                decIndent();
                appendIndent();
            }
        }
        endStatement(sid);
        appendTrailingComment(owner);
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
        var dangling = comments.dangling(node);
        var inlineDangling = !dangling.isEmpty() && dangling.stream().allMatch(c -> c.text.startsWith("/*") && c.text.indexOf('\n') < 0);
        var hasComments = !dangling.isEmpty() || comments.hasComments(code);
        if( code.getStatements().size() == 0 && (dangling.isEmpty() || inlineDangling) ) {
            for( var comment : dangling ) {
                comments.markEmitted(comment);
                append(' ');
                append(comment.text);
            }
            append(" }");
        }
        else if( code.getStatements().size() == 1 && code.getStatements().get(0) instanceof ExpressionStatement es && !shouldWrapExpression(node) && !hasComments && !comments.hasComments(es) ) {
            append(' ');
            visitStatementLabels(es);
            visit(es.getExpression());
            append(" }");
        }
        else {
            appendNewLine();
            incIndent();
            visit(code);
            appendDanglingComments(node);
            decIndent();
            appendIndent();
            append('}');
        }
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
        var wrap = hasTrailingComma(node) || hasComments(node, node.getExpressions());
        append('(');
        if( wrap )
            incIndent();
        visitPositionalArgs(node.getExpressions(), wrap);
        if( wrap ) {
            appendNewLine();
            appendDanglingComments(node);
            decIndent();
            appendIndent();
        }
        append(')');
    }

    @Override
    public void visitListExpression(ListExpression node) {
        var wrap = hasTrailingComma(node) || hasComments(node, node.getExpressions());
        append('[');
        if( wrap )
            incIndent();
        visitPositionalArgs(node.getExpressions(), wrap);
        if( wrap ) {
            appendNewLine();
            appendDanglingComments(node);
            decIndent();
            appendIndent();
        }
        append(']');
    }

    protected void visitPositionalArgs(List<Expression> args, boolean wrap) {
        var comma = wrap ? "," : ", ";
        var trailingComma = wrap && args.size() > 1;
        for( int i = 0; i < args.size(); i++ ) {
            var arg = args.get(i);
            if( wrap ) {
                appendNewLine();
                appendLeadingComments(arg);
                appendIndent();
            }
            visit(arg);
            if( trailingComma || i + 1 < args.size() )
                append(comma);
            if( wrap )
                appendTrailingComment(arg);
        }
    }

    @Override
    public void visitMapExpression(MapExpression node) {
        if( node.getMapEntryExpressions().isEmpty() ) {
            append("[:]");
            return;
        }
        var wrap = hasTrailingComma(node) || hasComments(node, node.getMapEntryExpressions());
        append('[');
        if( wrap )
            incIndent();
        visitNamedArgs(node.getMapEntryExpressions(), wrap);
        if( wrap ) {
            appendNewLine();
            appendDanglingComments(node);
            decIndent();
            appendIndent();
        }
        append(']');
    }

    protected void visitNamedArgs(List<MapEntryExpression> args, boolean wrap) {
        var comma = wrap ? "," : ", ";
        var trailingComma = wrap && args.size() > 1;
        for( int i = 0; i < args.size(); i++ ) {
            var arg = args.get(i);
            if( wrap ) {
                appendNewLine();
                appendLeadingComments(arg);
                appendIndent();
            }
            visit(arg);
            if( trailingComma || i + 1 < args.size() )
                append(comma);
            if( wrap )
                appendTrailingComment(arg);
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
            append(isTripleQuoted(text) ? reindentString(text) : text);
        else
            append(node.getText());
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
        var tripleQuoted = TDQ_STR.equals(quoteChar);
        append(quoteChar);
        var ss = node.getStrings();
        var vs = node.getValues();
        for( int i = 0; i < ss.size(); i++ ) {
            var string = ss.get(i);
            var text = (String) string.getNodeMetaData(ASTNodeMarker.VERBATIM_TEXT);
            if( text != null )
                append(tripleQuoted ? reindentString(text) : text);
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
        super.visit(node);
        if( number != null && number.intValue() > 0 )
            append(')');
    }

    // helpers

    private boolean hasComments(ASTNode container, List<? extends Expression> args) {
        if( !comments.dangling(container).isEmpty() )
            return true;
        for( var arg : args ) {
            if( comments.hasComments(arg) )
                return true;
        }
        return false;
    }

    /**
     * Determine whether any link of a method chain carries a comment.
     */
    private boolean hasChainComments(MethodCallExpression node) {
        Expression expr = node;
        while( expr instanceof MethodCallExpression mce && !mce.isImplicitThis() ) {
            var receiver = mce.getObjectExpression();
            if( !comments.leading(mce).isEmpty() || !comments.trailing(receiver).isEmpty() )
                return true;
            expr = receiver;
        }
        return false;
    }

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

    private static boolean isTripleQuoted(String text) {
        return text.startsWith(TDQ_STR) || text.startsWith(TSQ_STR);
    }

    private boolean shouldWrapExpression(Expression node) {
        return node.getLineNumber() < node.getLastLineNumber();
    }

    private boolean shouldWrapMethodCall(MethodCallExpression node) {
        if( hasTrailingComma(node.getArguments()) )
            return true;
        var start = node.getMethod();
        var end = node.getArguments();
        return start.getLineNumber() < end.getLastLineNumber();
    }

    private boolean shouldWrapMethodChain(MethodCallExpression node) {
        if( currentRootExpr != node )
            return false;
        // a comment between two links can only be emitted by wrapping
        if( hasChainComments(node) )
            return true;
        if( !shouldWrapExpression(node) )
            return false;

        Expression root = node;
        int depth = 0;
        while( root instanceof MethodCallExpression mce && !mce.isImplicitThis() ) {
            root = mce.getObjectExpression();
            depth += 1;
        }

        return shouldWrapExpression(root)
            ? false
            : depth >= 2;
    }

    private static final String SLASH_STR = "/";
    private static final String TDQ_STR = "\"\"\"";
    private static final String TSQ_STR = "'''";
    private static final String SQ_STR = "'";
    private static final String DQ_STR = "\"";

}
