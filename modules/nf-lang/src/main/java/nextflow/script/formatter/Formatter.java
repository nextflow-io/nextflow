/*
 * Copyright 2024-2025, Seqera Labs
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

import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.parser.CommentWriter;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.CodeVisitorSupport;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.expr.BinaryExpression;
import org.codehaus.groovy.ast.expr.BitwiseNegationExpression;
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

    private CommentWriter commentWriter;

    private StringBuilder builder = new StringBuilder();

    private int indentCount = 0;

    public Formatter(FormattingOptions options) {
        this.options = options;
        this.commentWriter = null;
    }

    public Formatter(FormattingOptions options, CommentWriter commentWriter) {
        this.options = options;
        this.commentWriter = commentWriter;
    }

    public void append(char c) {
        builder.append(c);
    }

    public void append(String str) {
        builder.append(str);
    }

    private char peekBuilder() {
        return builder.charAt(builder.length() - 1);
    }

    public void appendSpace() {
        if (Character.isWhitespace(peekBuilder())) {
            // We never want to append double space, or a space directly after an indent of a newline 
            return;
        } else {
            append(' ');
        }
    }

    public void appendPadding(int padding) {
        append(" ".repeat(padding));
    }

    public void appendIndent() {}

    public void _appendIndent() {
        var str = options.insertSpaces()
            ? " ".repeat(options.tabSize() * indentCount)
            : "\t".repeat(indentCount);
        builder.append(str);
    }

    boolean isFirst = true;
    private void appendNL() {
        if (isFirst) {
            isFirst = false;
            return;
        }
        builder.append('\n');
    }

    private void appendNLs(int i) {
        builder.append("\n".repeat(i));
    }

    public void appendNewLine() {
        appendNL();
        _appendIndent();
    }

    public void appendBlanks(int i) {
        appendNLs(i);
    }

    /*
     * Write a single comment while keeping track of breaking lines and WS before and after
     */
    public boolean writeComment(
        CommentWriter.Comment comment,
        boolean newlineBefore,
        boolean newlineAfter,
        boolean makeWSBefore,
        boolean makeWSAfter
    ) {
        var commentAndIsSLC = comment.write(); 
        String commentText = commentAndIsSLC.getV1();
        boolean isSLC = commentAndIsSLC.getV2();

        if( newlineBefore ) {
            appendNewLine();
        } else if( makeWSBefore ) {
            appendSpace();
        }

        append(commentText);

        if( isSLC || newlineAfter ) {
            return true;
        } 

        if( makeWSAfter ) {
            appendSpace();
        }

        return false;
    }

    public void writeComments(
        Map<String, List<CommentWriter.Comment>> comments,
        String key,
        boolean newlineBefore,
        boolean makeWSBefore,
        boolean makeWSAfter,
        boolean MLCOnOwnLines,
        boolean inExpression
    ) {
        if (comments.containsKey(key)) {
            writeComments(
                comments.get(key),
                newlineBefore,
                makeWSBefore,
                makeWSAfter,
                MLCOnOwnLines,
                inExpression
            );
        }
    }

    public void writeComments(
        List<CommentWriter.Comment> comments,
        boolean newlineBefore,
        boolean makeWSBefore,
        boolean makeWSAfter,
        boolean MLCOnOwnLines,
        boolean inExpression
    ) {
        if (comments.isEmpty()) return;
        /*
         * These are are the four cases:
         * - WSBefore && WSAfter:
         *   The first comment writes WS before and after,
         *   the rest are responsible for writing the WS after
         * 
         * - WSBefore = true, WSAfter = false:
         *   Each comment should only write WS before
         * 
         * - WSBefore = false, WSAfter = true:
         *   Each comment should only write WS after
         * 
         * 
         * - WSBefore = false, WSAfter = false:
         *   The first comment writes no whitespace at all (if it is not alone),
         *   the rest of the comment only write WS before
         */
        Iterator<CommentWriter.Comment> it = comments.iterator();
        var breakNext = writeComment(
            it.next(), newlineBefore, MLCOnOwnLines, makeWSBefore, makeWSAfter
        );
        while( it.hasNext() ) 
            breakNext = writeComment(
                it.next(), breakNext, MLCOnOwnLines, !makeWSAfter, makeWSAfter
            );

        if (breakNext && inExpression) {
            // If we inside an expression where the next token will not break the line
            // we are required to break the line here already
            appendNewLine(); 
        }
    }

    public void appendLeadingComments(
        ASTNode node,
        boolean newlineBefore,
        boolean makeWSBefore,
        boolean makeWSAfter,
        boolean MLCOnOwnLines,
        boolean inExpression
    ) {
        var comments = commentWriter.getLeadingComments(node);
        writeComments(
            comments,
            newlineBefore,
            makeWSBefore,
            makeWSAfter,
            MLCOnOwnLines,
            inExpression
        );
    }

    public void appendLeadingComments(ASTNode node) {
        appendLeadingComments(node, true, false, false, true, false);
    }

    public void appendLeadingCommentsExpression(ASTNode node) {
        appendLeadingComments(node, false, false, true, false, true);

    }

    public void appendLeadingCommentsStatement(ASTNode node) {
        appendLeadingComments(node, true, false, false, true, false);
    }

    public void appendTrailingComments(
        ASTNode node,
        boolean newlineBefore,
        boolean makeWSBefore,
        boolean makeWSAfter,
        boolean MLCOnOwnLines,
        boolean inExpression
    ) {
        var comments = commentWriter.getTrailingComments(node);
        if (comments.containsKey("STANDARD")) {
            writeComments(
                comments.get("STANDARD"),
                newlineBefore,
                makeWSBefore,
                makeWSAfter,
                MLCOnOwnLines,
                inExpression               
            );
        }
    }

    public void appendTrailingInside(
        Map<String, List<CommentWriter.Comment>> comments,
        String key
    ) {
        writeComments(comments, key, false, true, true, false, true);
    }

    public void appendWithinComments(
        Map<String, List<CommentWriter.Comment>> comments,
        String key
    ) {
        writeComments(comments, key, false, true, false, false, true);
    }

    public void appendTrailingComments(ASTNode node) {
        appendTrailingComments(node, false, true, false, false, false);
    }

    public void appendTrailingCommentsExpression(ASTNode node) {
        appendTrailingComments(node, false, true, true, false, true);
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
    public void visitIfElse(IfStatement node) {
        visitIfElse(node, true);
    }

    protected void visitIfElse(IfStatement node, boolean preIndent) {
        appendLeadingComments(node);
        appendNewLine();
        append("if");
        appendSpace();
        append('(');
        visit(node.getBooleanExpression());
        append(')');
        appendSpace();
        append("{");
        incIndent();
        visit(node.getIfBlock());
        decIndent();
        appendNewLine();
        append("}");
        if( node.getElseBlock() instanceof IfStatement is ) {
            appendNewLine();
            append("else");
            appendSpace();
            visitIfElse(is, false);
        }
        else if( !(node.getElseBlock() instanceof EmptyStatement) ) {
            appendNewLine();
            append("else");
            appendSpace();
            append('{');
            incIndent();
            visit(node.getElseBlock());
            decIndent();
            appendNewLine();
            append("}");
        }
    }

    private Expression currentRootExpr;

    @Override
    public void visitExpressionStatement(ExpressionStatement node) {
        var cre = currentRootExpr;
        currentRootExpr = node.getExpression();
        appendLeadingComments(node);
        appendNewLine();
        visitStatementLabels(node);
        visit(node.getExpression());
        appendTrailingComments(node);
        currentRootExpr = cre;
    }

    private void visitStatementLabels(ExpressionStatement node) {
        if( node.getStatementLabels() == null )
            return;
        for( var label : DefaultGroovyMethods.asReversed(node.getStatementLabels()) ) {
            append(label);
            append(":");
            appendSpace();
        }
    }

    @Override
    public void visitReturnStatement(ReturnStatement node) {
        var cre = currentRootExpr;
        currentRootExpr = node.getExpression();
        appendLeadingComments(node);
        appendNewLine();
        append("return");
        appendSpace();
        visit(node.getExpression());
        currentRootExpr = cre;
    }

    @Override
    public void visitAssertStatement(AssertStatement node) {
        appendLeadingComments(node);
        appendNewLine();
        append("assert");
        appendSpace();
        visit(node.getBooleanExpression());
        if( !(node.getMessageExpression() instanceof ConstantExpression ce && ce.isNullExpression()) ) {
            appendSpace();
            append(":");
            appendSpace();
            visit(node.getMessageExpression());
        }
    }

    @Override
    public void visitTryCatchFinally(TryCatchStatement node) {
        appendLeadingComments(node);
        appendNewLine();
        append("try");
        appendSpace();
        append('{');
        incIndent();
        appendNewLine();
        visit(node.getTryStatement());
        decIndent();
        appendNewLine();
        append("}");
        for( var catchStatement : node.getCatchStatements() ) {
            visit(catchStatement);
        }
    }

    @Override
    public void visitThrowStatement(ThrowStatement node) {
        appendLeadingComments(node);
        appendNewLine();
        append("throw");
        appendSpace();
        visit(node.getExpression());
    }

    @Override
    public void visitCatchStatement(CatchStatement node) {
        appendLeadingComments(node);
        appendNewLine();
        append("catch");
        appendSpace();
        append('(');

        var variable = node.getVariable();
        var type = variable.getType();
        if( !ClassHelper.isObjectType(type) ) {
            append(type.getNameWithoutPackage());
            appendSpace();
        }
        append(variable.getName());
        append(')');
        appendSpace();
        append('{');
        incIndent();
        visit(node.getCode());
        decIndent();
        appendNewLine();
        append("}");
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
                if( !nextflow.script.types.Types.isNamespace(receiver.getType()) ) {
                    appendNewLine();
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
                decIndent();
                appendNewLine();
            }
            append(')');
        }
        if( lastClosureArg ) {
            appendSpace();
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
        append("new");
        appendSpace();
        visitTypeAnnotation(node.getType());
        append('(');
        visitArguments(asMethodCallArguments(node), false);
        append(')');
    }

    public void visitDirective(MethodCallExpression call) {
        appendLeadingComments(call);
        appendNewLine();
        append(call.getMethodAsString());
        var arguments = asMethodCallArguments(call);
        if( !arguments.isEmpty() ) {
            appendSpace();
            visitArguments(arguments, false);
        }
    }

    public void visitArguments(List<Expression> args, boolean wrap) {
        var hasNamedArgs = args.size() > 0 && args.get(0) instanceof NamedArgumentListExpression;
        var positionalArgs = hasNamedArgs
            ? DefaultGroovyMethods.tail(args)
            : args;
        visitPositionalArgs(positionalArgs, wrap);
        if( hasNamedArgs ) {
            if( positionalArgs.size() > 0 )
                append(',');
                if (!wrap) appendSpace();
            var mapX = (MapExpression)args.get(0);
            visitNamedArgs(mapX.getMapEntryExpressions(), wrap);
        }
    }

    private boolean inVariableDeclaration;

    private boolean inWrappedPipeChain;

    @Override
    public void visitBinaryExpression(BinaryExpression node) {
        appendLeadingCommentsExpression(node);
        if( node instanceof DeclarationExpression ) {
            appendNewLine(); // A declaration expression is treated as a statement here, so break the line
            append("def");
            appendSpace();
            inVariableDeclaration = true;
            visit(node.getLeftExpression());
            inVariableDeclaration = false;
            var source = node.getRightExpression();
            if( !(source instanceof EmptyExpression) ) {
                append("=");
                appendSpace();
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

        var beginWrappedPipeChain = shouldWrapPipeExpression(node);
        if( beginWrappedPipeChain )
            inWrappedPipeChain = true;

        visit(node.getLeftExpression());

        if( inWrappedPipeChain ) {
            incIndent();
            appendNewLine();
        }
        else {
            appendSpace();
        }
        append(node.getOperation().getText());
        appendSpace();

        var iwpc = inWrappedPipeChain;
        inWrappedPipeChain = false;
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
        inWrappedPipeChain = iwpc;

        if( inWrappedPipeChain )
            decIndent();

        if( beginWrappedPipeChain )
            inWrappedPipeChain = false;
    }

    @Override
    public void visitTernaryExpression(TernaryExpression node) {
        appendLeadingCommentsExpression(node);
        if( shouldWrapExpression(node) ) {
            visit(node.getBooleanExpression());
            incIndent();
            appendNewLine();
            append("?");
            appendSpace();
            visit(node.getTrueExpression());
            appendNewLine();
            append(":");
            appendSpace();
            visit(node.getFalseExpression());
            decIndent();
        }
        else {
            visit(node.getBooleanExpression());
            appendSpace();
            append("?");
            appendSpace();
            visit(node.getTrueExpression());
            appendSpace();
            append(":");
            appendSpace();
            visit(node.getFalseExpression());
        }
    }

    @Override
    public void visitShortTernaryExpression(ElvisOperatorExpression node) {
        appendLeadingCommentsExpression(node);
        visit(node.getTrueExpression());
        appendSpace();
        append("?:");
        appendSpace();
        visit(node.getFalseExpression());
    }

    @Override
    public void visitNotExpression(NotExpression node) {
        appendLeadingCommentsExpression(node);
        append('!');
        visit(node.getExpression());
    }

    @Override
    public void visitClosureExpression(ClosureExpression node) {
        appendLeadingCommentsExpression(node);
        append('{');
        if( node.getParameters() != null && node.getParameters().length > 0 ) {
            appendSpace();
            visitParameters(node.getParameters());
            appendSpace();
            append("->");
        }
        var code = (BlockStatement) node.getCode();
        if( code.getStatements().size() == 0 ) {
            appendSpace();
            append("}");
        }
        else if( code.getStatements().size() == 1 && code.getStatements().get(0) instanceof ExpressionStatement es && !shouldWrapExpression(node) ) {
            appendSpace();
            visitStatementLabels(es);
            visit(es.getExpression());
            appendSpace();
            append("}");
        }
        else {
            incIndent();
            visit(code);
            decIndent();
            appendNewLine();
            append('}');
        }
    }

    public void visitParameters(Parameter[] parameters) {
        for( int i = 0; i < parameters.length; i++ ) {
            var param = parameters[i];
            if( isLegacyType(param.getType()) ) {
                visitTypeAnnotation(param.getType());
                appendSpace();
            }
            append(param.getName());
            if( param.hasInitialExpression() ) {
                appendSpace();
                append("=");
                appendSpace();
                visit(param.getInitialExpression());
            }
            if( i + 1 < parameters.length ) {
                append(",");
                appendSpace();
            }
        }
    }

    @Override
    public void visitTupleExpression(TupleExpression node) {
        appendLeadingCommentsExpression(node);
        var wrap = shouldWrapExpression(node);
        append('(');
        if( wrap )
            incIndent();
        visitPositionalArgs(node.getExpressions(), wrap);
        if( wrap ) {
            decIndent();
            appendNewLine();
        }
        append(')');
    }

    @Override
    public void visitListExpression(ListExpression node) {
        appendLeadingCommentsExpression(node);
        var wrap = hasTrailingComma(node) || shouldWrapExpression(node);
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
            }
            visit(args.get(i));
            if( trailingComma || i + 1 < args.size() )
                append(comma);
        }
    }

    @Override
    public void visitMapExpression(MapExpression node) {
        appendLeadingCommentsExpression(node);
        if( node.getMapEntryExpressions().isEmpty() ) {
            append("[:]");
            return;
        }
        var wrap = hasTrailingComma(node) || shouldWrapExpression(node);
        append('[');
        if( wrap )
            incIndent();
        visitNamedArgs(node.getMapEntryExpressions(), wrap);
        if( wrap ) {
            decIndent();
            appendNewLine();
        }
        append(']');
    }

    protected void visitNamedArgs(List<MapEntryExpression> args, boolean wrap) {
        var comma = wrap ? "," : ", ";
        var trailingComma = wrap && args.size() > 1;
        for( int i = 0; i < args.size(); i++ ) {
            if( wrap ) {
                appendNewLine();
            }
            visit(args.get(i));
            if( trailingComma || i + 1 < args.size() )
                append(comma);
        }
    }

    @Override
    public void visitMapEntryExpression(MapEntryExpression node) {
        appendLeadingCommentsExpression(node);
        visit(node.getKeyExpression());
        append(":");
        appendSpace();
        visit(node.getValueExpression());
    }

    @Override
    public void visitRangeExpression(RangeExpression node) {
        appendLeadingCommentsExpression(node);
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
        appendLeadingCommentsExpression(node);
        append('-');
        visit(node.getExpression());
    }

    @Override
    public void visitUnaryPlusExpression(UnaryPlusExpression node) {
        appendLeadingCommentsExpression(node);
        append('+');
        visit(node.getExpression());
    }

    @Override
    public void visitBitwiseNegationExpression(BitwiseNegationExpression node) {
        appendLeadingCommentsExpression(node);
        append('~');
        visit(node.getExpression());
    }

    @Override
    public void visitCastExpression(CastExpression node) {
        appendLeadingCommentsExpression(node);
        visit(node.getExpression());
        appendSpace();
        append("as");
        appendSpace();
        visitTypeAnnotation(node.getType());
    }

    @Override
    public void visitConstantExpression(ConstantExpression node) {
        appendLeadingCommentsExpression(node);
        var text = (String) node.getNodeMetaData(ASTNodeMarker.VERBATIM_TEXT);
        if( text != null )
            append(text);
        else
            append(node.getText());
        appendTrailingCommentsExpression(node);
    }

    @Override
    public void visitClassExpression(ClassExpression node) {
        appendLeadingCommentsExpression(node);
        visitTypeAnnotation(node.getType());
    }

    public void visitTypeAnnotation(ClassNode type) {
        appendLeadingCommentsExpression(type);
        if( isLegacyType(type) ) {
            append(type.getNodeMetaData(ASTNodeMarker.LEGACY_TYPE));
            return;
        }

        append(nextflow.script.types.Types.getName(type));
    }

    @Override
    public void visitVariableExpression(VariableExpression node) {
        appendLeadingCommentsExpression(node);
        if( inVariableDeclaration && isLegacyType(node.getType()) ) {
            visitTypeAnnotation(node.getType());
            appendSpace();
        }
        append(node.getText());
        appendTrailingCommentsExpression(node);
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
        appendLeadingCommentsExpression(node);
        // see also: GStringUtil.writeToImpl()
        var quoteChar = (String) node.getNodeMetaData(ASTNodeMarker.QUOTE_CHAR, k -> DQ_STR);
        append(quoteChar);
        var ss = node.getStrings();
        var vs = node.getValues();
        for( int i = 0; i < ss.size(); i++ ) {
            var string = ss.get(i);
            if( string.getNodeMetaData(ASTNodeMarker.VERBATIM_TEXT) != null )
                visit(string);
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

    private static boolean hasTrailingComma(Expression node) {
        return node.getNodeMetaData(ASTNodeMarker.TRAILING_COMMA) != null;
    }

    public static boolean isLegacyType(ClassNode cn) {
        return cn.getNodeMetaData(ASTNodeMarker.LEGACY_TYPE) != null;
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

    private boolean shouldWrapPipeExpression(BinaryExpression node) {
        return currentRootExpr == node && node.getOperation().isA(Types.PIPE) && shouldWrapExpression(node);
    }

    private static final String SLASH_STR = "/";
    private static final String TDQ_STR = "\"\"\"";
    private static final String TSQ_STR = "'''";
    private static final String SQ_STR = "'";
    private static final String DQ_STR = "\"";

}
