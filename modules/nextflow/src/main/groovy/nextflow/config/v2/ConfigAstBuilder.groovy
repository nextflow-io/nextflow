/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.config.v2

import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.antlr.ConfigLexer
import nextflow.antlr.ConfigParser
import nextflow.antlr.DescriptiveErrorStrategy
import org.antlr.v4.runtime.ANTLRErrorListener
import org.antlr.v4.runtime.CharStream
import org.antlr.v4.runtime.CharStreams
import org.antlr.v4.runtime.CommonTokenStream
import org.antlr.v4.runtime.Parser
import org.antlr.v4.runtime.ParserRuleContext
import org.antlr.v4.runtime.RecognitionException
import org.antlr.v4.runtime.Recognizer
import org.antlr.v4.runtime.Token as ParserToken
import org.antlr.v4.runtime.atn.ATNConfigSet
import org.antlr.v4.runtime.atn.PredictionMode
import org.antlr.v4.runtime.dfa.DFA
import org.antlr.v4.runtime.misc.ParseCancellationException
import org.antlr.v4.runtime.tree.TerminalNode
import org.apache.groovy.parser.antlr4.GroovySyntaxError
import org.apache.groovy.parser.antlr4.util.StringUtils
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassHelper
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.GenericsType
import org.codehaus.groovy.ast.ModuleNode
import org.codehaus.groovy.ast.NodeMetaDataHandler
import org.codehaus.groovy.ast.Parameter
import org.codehaus.groovy.ast.VariableScope
import org.codehaus.groovy.ast.expr.ArgumentListExpression
import org.codehaus.groovy.ast.expr.BinaryExpression
import org.codehaus.groovy.ast.expr.BitwiseNegationExpression
import org.codehaus.groovy.ast.expr.BooleanExpression
import org.codehaus.groovy.ast.expr.ClassExpression
import org.codehaus.groovy.ast.expr.ConstantExpression
import org.codehaus.groovy.ast.expr.EmptyExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.GStringExpression
import org.codehaus.groovy.ast.expr.MapExpression
import org.codehaus.groovy.ast.expr.MapEntryExpression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.NotExpression
import org.codehaus.groovy.ast.expr.PostfixExpression
import org.codehaus.groovy.ast.expr.PrefixExpression
import org.codehaus.groovy.ast.expr.PropertyExpression
import org.codehaus.groovy.ast.expr.RangeExpression
import org.codehaus.groovy.ast.expr.SpreadExpression
import org.codehaus.groovy.ast.expr.SpreadMapExpression
import org.codehaus.groovy.ast.expr.TupleExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.ast.expr.UnaryMinusExpression
import org.codehaus.groovy.ast.expr.UnaryPlusExpression
import org.codehaus.groovy.ast.stmt.AssertStatement
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.EmptyStatement
import org.codehaus.groovy.ast.stmt.ReturnStatement
import org.codehaus.groovy.ast.stmt.Statement
import org.codehaus.groovy.control.CompilationFailedException
import org.codehaus.groovy.control.CompilePhase
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.control.messages.SyntaxErrorMessage
import org.codehaus.groovy.syntax.Numbers
import org.codehaus.groovy.syntax.SyntaxException
import org.codehaus.groovy.syntax.Token
import org.codehaus.groovy.syntax.Types

import static nextflow.antlr.ConfigParser.*
import static nextflow.antlr.PositionConfigureUtils.configureAST as ast
import static nextflow.ast.ASTHelpers.*
import static org.codehaus.groovy.ast.tools.GeneralUtils.*

@Slf4j
@CompileStatic
class ConfigAstBuilder {

    private SourceUnit sourceUnit
    private ModuleNode moduleNode
    private ConfigLexer lexer
    private ConfigParser parser

    private Tuple2<ParserRuleContext,Exception> numberFormatError

    ConfigAstBuilder(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit
        this.moduleNode = new ModuleNode(sourceUnit)

        final charStream = createCharStream(sourceUnit)
        this.lexer = new ConfigLexer(charStream)
        this.parser = new ConfigParser(new CommonTokenStream(lexer))
        parser.setErrorHandler(new DescriptiveErrorStrategy(charStream))
    }

    private CharStream createCharStream(SourceUnit sourceUnit) {
        try {
            return CharStreams.fromReader(
                    new BufferedReader(sourceUnit.getSource().getReader()),
                    sourceUnit.getName())
        } catch (IOException e) {
            throw new RuntimeException("Error occurred when reading source code.", e)
        }
    }

    private CompilationUnitContext buildCST() {
        try {
            final tokenStream = parser.getInputStream()
            try {
                return buildCST(PredictionMode.SLL)
            }
            catch( Throwable t ) {
                // if some syntax error occurred in the lexer, no need to retry the powerful LL mode
                if( t instanceof GroovySyntaxError && t.getSource() == GroovySyntaxError.LEXER )
                    throw t

                log.trace "Parsing mode SLL failed, falling back to LL"
                tokenStream.seek(0)
                return buildCST(PredictionMode.LL)
            }
        }
        catch( Throwable t ) {
            throw convertException(t)
        }
    }

    private CompilationUnitContext buildCST(PredictionMode predictionMode) {
        parser.getInterpreter().setPredictionMode(predictionMode)

        if( predictionMode == PredictionMode.SLL )
            removeErrorListeners()
        else
            addErrorListeners()

        return parser.compilationUnit()
    }

    private CompilationFailedException convertException(Throwable t) {
        if( t instanceof CompilationFailedException )
            return t
        else if( t instanceof ParseCancellationException )
            return createParsingFailedException(t.getCause())
        else
            return createParsingFailedException(t)
    }

    ModuleNode buildAST(SourceUnit sourceUnit) {
        try {
            return compilationUnit(buildCST())
        } catch (Throwable t) {
            throw convertException(t)
        }
    }

    /// CONFIG STATEMENTS

    private ModuleNode compilationUnit(CompilationUnitContext ctx) {
        for( final stmt : ctx.configStatement() )
            moduleNode.addStatement(configStatement(stmt))

        if( moduleNode.isEmpty() )
            moduleNode.addStatement(ReturnStatement.RETURN_NULL_OR_VOID)

        final scriptClassNode = moduleNode.getScriptClassDummy()
        final statements = moduleNode.getStatementBlock().getStatements()
        if( scriptClassNode && !statements.isEmpty() ) {
            final first = statements.first()
            final last = statements.last()
            scriptClassNode.setSourcePosition(first)
            scriptClassNode.setLastColumnNumber(last.getLastColumnNumber())
            scriptClassNode.setLastLineNumber(last.getLastLineNumber())
        }

        if( numberFormatError != null )
            throw createParsingFailedException(numberFormatError.getV2().getMessage(), numberFormatError.getV1())

        return moduleNode
    }

    private Statement configStatement(ConfigStatementContext ctx) {
        if( ctx instanceof ConfigIncludeStmtAltContext )
            return ast( configInclude(ctx.configInclude()), ctx )

        if( ctx instanceof ConfigAssignmentStmtAltContext )
            return ast( configAssignment(ctx.configAssignment()), ctx )

        if( ctx instanceof ConfigBlockStmtAltContext )
            return ast( configBlock(ctx.configBlock()), ctx )

        throw createParsingFailedException("Invalid config statement: ${ctx.text}", ctx)
    }

    private Statement configInclude(ConfigIncludeContext ctx) {
        final source = expression(ctx.expression())
        stmt(ast( callThisX('includeConfig', args(source)), ctx ))
    }

    private Statement configAssignment(ConfigAssignmentContext ctx) {
        final names = listX( ctx.configPathExpression().identifier().collect( this.&configPathElement ) )
        final right = expression(ctx.expression())
        stmt(ast( callThisX('assign', args(names, right)), ctx ))
    }

    private Expression configPathElement(IdentifierContext ctx) {
        ast( constX(identifier(ctx)), ctx )
    }

    private Statement configBlock(ConfigBlockContext ctx) {
        final name = ctx.identifier()
            ? constX(identifier(ctx.identifier()))
            : constX(stringLiteral(ctx.stringLiteral()))
        final statements = ctx.configBlockStatement().collect( this.&configBlockStatement )
        final closure = closureX(new BlockStatement(statements, new VariableScope()))
        stmt(ast( callThisX('block', args(name, closure)), ctx ))
    }

    private Statement configBlockStatement(ConfigBlockStatementContext ctx) {
        if( ctx instanceof ConfigIncludeBlockStmtAltContext )
            return ast( configInclude(ctx.configInclude()), ctx )

        if( ctx instanceof ConfigAssignmentBlockStmtAltContext )
            return ast( configAssignment(ctx.configAssignment()), ctx )

        if( ctx instanceof ConfigBlockBlockStmtAltContext )
            return ast( configBlock(ctx.configBlock()), ctx )

        if( ctx instanceof ConfigSelectorBlockStmtAltContext )
            return ast( configSelector(ctx.configSelector()), ctx )

        throw createParsingFailedException("Invalid statement in config block: ${ctx.text}", ctx)
    }

    private Statement configSelector(ConfigSelectorContext ctx) {
        final kind = ctx.kind.text
        final target = configSelectorTarget(ctx.target)
        final statements = ctx.configAssignment().collect( this.&configAssignment )
        final closure = closureX(new BlockStatement(statements, new VariableScope()))
        stmt(ast( callThisX(kind, args(target, closure)), ctx ))
    }

    private Expression configSelectorTarget(ConfigSelectorTargetContext ctx) {
        ctx.identifier()
            ? constX(identifier(ctx.identifier()))
            : constX(stringLiteral(ctx.stringLiteral()))
    }

    /// GROOVY STATEMENTS

    private Statement statement(StatementContext ctx) {
        if( ctx instanceof ReturnStmtAltContext )
            return ast( returnStatement(ctx.expression()), ctx )

        if( ctx instanceof AssertStmtAltContext )
            return ast( assertStatement(ctx.assertStatement()), ctx )

        if( ctx instanceof VariableDeclarationStmtAltContext )
            return ast( variableDeclaration(ctx.variableDeclaration()), ctx )

        if( ctx instanceof MultipleAssignmentStmtAltContext )
            return ast( assignment(ctx.multipleAssignmentStatement()), ctx )

        if( ctx instanceof AssignmentStmtAltContext )
            return ast( assignment(ctx.assignmentStatement()), ctx )

        if( ctx instanceof ExpressionStmtAltContext )
            return ast( expressionStatement(ctx.expressionStatement()), ctx )

        if( ctx instanceof EmptyStmtAltContext )
            return new EmptyStatement()

        throw createParsingFailedException("Invalid Groovy statement: ${ctx.text}", ctx)
    }

    private Statement returnStatement(ExpressionContext ctx) {
        final result = ctx
            ? expression(ctx)
            : ConstantExpression.EMPTY_EXPRESSION
        returnS(result)
    }

    private Statement assertStatement(AssertStatementContext ctx) {
        final condition = ast( new BooleanExpression(expression(ctx.condition)), ctx.condition )
        ctx.message
            ? new AssertStatement(condition, expression(ctx.message))
            : new AssertStatement(condition)
    }

    private Statement variableDeclaration(VariableDeclarationContext ctx) {
        if( ctx.typeNamePairs() ) {
            // multiple assignment
            final variables = ctx.typeNamePairs().typeNamePair().collect { pair ->
                final name = identifier(pair.identifier())
                final type = type(pair.type())
                ast( varX(name, type), pair )
            }
            final target = variables.size() > 1
                ? new ArgumentListExpression(variables as List<Expression>)
                : variables.first()
            final initializer = expression(ctx.initializer)
            return stmt(ast( declX(target, initializer), ctx ))
        }
        else {
            // single assignment
            final type = type(ctx.type())
            final decl = ctx.variableDeclarator()
            final name = identifier(decl.identifier())
            final target = ast( varX(name, type), ctx )
            final initializer = decl.initializer
                ? expression(decl.initializer)
                : EmptyExpression.INSTANCE
            return stmt(ast( declX(target, initializer), ctx ))
        }
    }

    private Statement assignment(MultipleAssignmentStatementContext ctx) {
        final vars = ctx.variableNames().identifier().collect( this.&variableName )
        final left = ast( new TupleExpression(vars), ctx.variableNames() )
        final right = expression(ctx.right)
        return stmt(assignX(left, right))
    }

    private Expression variableName(IdentifierContext ctx) {
        return ast( varX(identifier(ctx)), ctx )
    }

    private Statement assignment(AssignmentStatementContext ctx) {
        final left = expression(ctx.left)
        if( left instanceof VariableExpression && isInsideParentheses(left) ) {
            if( left.<Number>getNodeMetaData(INSIDE_PARENTHESES_LEVEL).intValue() > 1 )
                throw createParsingFailedException("Nested parenthesis is not allowed in multiple assignment, e.g. ((a)) = b", ctx)

            final tuple = ast( new TupleExpression(left), ctx.left )
            return stmt(assignX(tuple, expression(ctx.right)))
        }

        if ( isAssignmentLhsValid(left) )
            return stmt(assignX(left, expression(ctx.right)))

        throw createParsingFailedException("The left-hand side of an assignment should be a variable or a property expression", ctx)
    }

    private boolean isAssignmentLhsValid(Expression left) {
        // e.g. p = 123
        if( left instanceof VariableExpression && !isInsideParentheses(left) )
            return true
        // e.g. obj.p = 123
        if( left instanceof PropertyExpression )
            return true
        // e.g. map[a] = 123 OR map['a'] = 123 OR map["$a"] = 123
        if( left instanceof BinaryExpression && left.operation.type == Types.LEFT_SQUARE_BRACKET )
            return true
        return false
    }

    private Statement expressionStatement(ExpressionStatementContext ctx) {
        final base = expression(ctx.expression())
        final expression = ctx.argumentList()
            ? methodCall(base, argumentList(ctx.argumentList()))
            : base
        return ast( stmt(expression), ctx )
    }

    /// GROOVY EXPRESSIONS

    private Expression expression(ExpressionContext ctx) {
        if( ctx instanceof AddSubExprAltContext )
            return ast( binary(ctx.left, ctx.op, ctx.right), ctx )

        if( ctx instanceof BitwiseAndExprAltContext )
            return ast( binary(ctx.left, ctx.op, ctx.right), ctx )

        if( ctx instanceof BitwiseOrExprAltContext )
            return ast( binary(ctx.left, ctx.op, ctx.right), ctx )

        if( ctx instanceof CastExprAltContext ) {
            final type = type(ctx.type())
            final operand = castOperand(ctx.castOperandExpression())
            return ast( castX(type, operand), ctx )
        }

        if( ctx instanceof ConditionalExprAltContext )
            return ast( ternary(ctx), ctx )

        if( ctx instanceof EqualityExprAltContext )
            return ast( binary(ctx.left, ctx.op, ctx.right), ctx )

        if( ctx instanceof ExclusiveOrExprAltContext )
            return ast( binary(ctx.left, ctx.op, ctx.right), ctx )

        if( ctx instanceof LogicalAndExprAltContext )
            return ast( binary(ctx.left, ctx.op, ctx.right), ctx )

        if( ctx instanceof LogicalOrExprAltContext )
            return ast( binary(ctx.left, ctx.op, ctx.right), ctx )

        if( ctx instanceof MultDivExprAltContext )
            return ast( binary(ctx.left, ctx.op, ctx.right), ctx )

        if( ctx instanceof PathExprAltContext )
            return path(ctx.pathExpression())

        if( ctx instanceof PostfixExprAltContext )
            return ast( postfix(ctx.pathExpression(), ctx.op), ctx )

        if( ctx instanceof PowerExprAltContext )
            return ast( binary(ctx.left, ctx.op, ctx.right), ctx )

        if( ctx instanceof PrefixExprAltContext )
            return ast( prefix(expression(ctx.expression()), ctx.op), ctx )

        if( ctx instanceof RegexExprAltContext )
            return ast( binary(ctx.left, ctx.op, ctx.right), ctx )

        if( ctx instanceof RelationalExprAltContext )
            return ast( binary(ctx.left, ctx.op, ctx.right), ctx )

        if( ctx instanceof RelationalCastExprAltContext ) {
            final operand = expression(ctx.expression())
            final type = type(ctx.type())
            return ast( asX(type, operand), ctx )
        }

        if( ctx instanceof RelationalTypeExprAltContext ) {
            final right = ast( new ClassExpression(type(ctx.type(), false)), ctx.type() )
            return ast( binary(ctx.left, ctx.op, right), ctx )
        }

        if( ctx instanceof ShiftExprAltContext )
            return ast( shift(ctx), ctx )

        if( ctx instanceof UnaryAddExprAltContext )
            return ast( unaryAdd(expression(ctx.expression()), ctx.op), ctx )

        if( ctx instanceof UnaryNotExprAltContext )
            return ast( unaryNot(expression(ctx.expression()), ctx.op), ctx )

        throw createParsingFailedException("Invalid Groovy expression: ${ctx.text}", ctx)
    }

    private Expression binary(ExpressionContext left, ParserToken op, ExpressionContext right) {
        binX(expression(left), token(op), expression(right))
    }

    private Expression binary(ExpressionContext left, ParserToken op, Expression right) {
        binX(expression(left), token(op), right)
    }

    private Expression castOperand(CastOperandExpressionContext ctx) {
        if( ctx instanceof CastCastExprAltContext ) {
            final type = type(ctx.type())
            final operand = castOperand(ctx.castOperandExpression())
            return ast( castX(type, operand), ctx )
        }

        if( ctx instanceof PathCastExprAltContext )
            return path(ctx.pathExpression())

        if( ctx instanceof PostfixCastExprAltContext )
            return ast( postfix(ctx.pathExpression(), ctx.op), ctx )

        if( ctx instanceof PrefixCastExprAltContext )
            return ast( prefix(castOperand(ctx.castOperandExpression()), ctx.op), ctx )

        if( ctx instanceof UnaryAddCastExprAltContext )
            return ast( unaryAdd(castOperand(ctx.castOperandExpression()), ctx.op), ctx )

        if( ctx instanceof UnaryNotCastExprAltContext )
            return ast( unaryNot(castOperand(ctx.castOperandExpression()), ctx.op), ctx )

        throw createParsingFailedException("Invalid Groovy expression: ${ctx.text}", ctx)
    }

    private Expression postfix(PathExpressionContext ctx, ParserToken op) {
        new PostfixExpression(path(ctx), token(op))
    }

    private Expression prefix(Expression expression, ParserToken op) {
        new PrefixExpression(token(op), expression)
    }

    private Expression shift(ShiftExprAltContext ctx) {
        final left = expression(ctx.left)
        final right = expression(ctx.right)

        if( ctx.riOp )
            return new RangeExpression(left, right, true)
        if( ctx.reOp )
            return new RangeExpression(left, right, false, true)

        def op
        if( ctx.dlOp )
            op = token(ctx.dlOp, 2)
        if( ctx.dgOp )
            op = token(ctx.dgOp, 2)
        if( ctx.tgOp )
            op = token(ctx.tgOp, 3)

        return binX(left, op, right)
    }

    private Expression ternary(ConditionalExprAltContext ctx) {
        if( ctx.ELVIS() )
            return elvisX(expression(ctx.condition), expression(ctx.fb))

        final condition = ast( boolX(expression(ctx.condition)), ctx.condition )
        return ternaryX(condition, expression(ctx.tb), expression(ctx.fb))
    }

    private Expression unaryAdd(Expression expression, ParserToken op) {
        if( op.type == ConfigParser.ADD )
            return new UnaryPlusExpression(expression)

        if( op.type == ConfigParser.SUB )
            return new UnaryMinusExpression(expression)

        throw new IllegalStateException()
    }

    private Expression unaryNot(Expression expression, ParserToken op) {
        if( op.type == ConfigParser.NOT )
            return new NotExpression(expression)

        if( op.type == ConfigParser.BITNOT )
            return new BitwiseNegationExpression(expression)

        throw new IllegalStateException()
    }

    /// -- PATH EXPRESSIONS

    private Expression path(PathExpressionContext ctx) {
        try {
            final primary = primary(ctx.primary())
            return ctx.pathElement().inject(primary, (acc, el) -> pathElement(acc, el))
        }
        catch( IllegalStateException e ) {
            throw createParsingFailedException("Invalid Groovy expression: ${ctx.text}", ctx)
        }
    }

    private Expression pathElement(Expression expression, PathElementContext ctx) {
        if( ctx instanceof PropertyPathExprAltContext )
            return ast( pathPropertyElement(expression, ctx), ctx )

        if( ctx instanceof ClosurePathExprAltContext )
            return ast( pathClosureElement(expression, ctx.closure()), ctx )

        if( ctx instanceof ArgumentsPathExprAltContext )
            return ast( pathArgumentsElement(expression, ctx.arguments().argumentList()), ctx )

        if( ctx instanceof IndexPathExprAltContext )
            return ast( pathIndexElement(expression, ctx.indexPropertyArgs()), ctx )

        throw new IllegalStateException()
    }

    private Expression pathPropertyElement(Expression expression, PropertyPathExprAltContext ctx) {
        final property = ast( constX(namePart(ctx.namePart())), ctx.namePart() )
        final safe = ctx.SAFE_DOT() != null || ctx.SPREAD_DOT() != null
        final result = new PropertyExpression(expression, property, safe)
        if( ctx.SPREAD_DOT() )
            result.setSpreadSafe(true)
        return result
    }

    private String namePart(NamePartContext ctx) {
        if( ctx.keywords() )
            return keywords(ctx.keywords())

        if( ctx.identifier() )
            return identifier(ctx.identifier())

        if( ctx.stringLiteral() )
            return stringLiteral(ctx.stringLiteral())

        throw new IllegalStateException()
    }

    private Expression pathClosureElement(Expression expression, ClosureContext ctx) {
        final closure = ast( closure(ctx), ctx )

        if( expression instanceof MethodCallExpression ) {
            // append closure to method call arguments
            final call = (MethodCallExpression)expression

            // normal arguments, e.g. 1, 2
            if ( call.arguments !instanceof ArgumentListExpression )
                throw new IllegalStateException()

            final arguments = (ArgumentListExpression)call.arguments
            arguments.addExpression(closure)
            return call

            // TODO: only needed if namedArgs uses TupleExpression
            // named arguments, e.g. x: 1, y: 2
            // if ( arguments instanceof TupleExpression ) {
            //     final tuple = (TupleExpression) arguments
            //     if( !tuple.expressions )
            //         throw new IllegalStateException()
            //     final namedArguments = (NamedArgumentListExpression) tuple.getExpression(0)
            //     call.arguments = args( mapX(namedArguments.mapEntryExpressions), closure )
            //     return call
            // }
        }

        final arguments = ast( args(closure), closure )

        // e.g. obj.m { }
        if( expression instanceof PropertyExpression )
            return propMethodCall(expression, arguments)

        // e.g. m { }, "$m" { }, "m" { }
        if( expression instanceof VariableExpression || expression instanceof GStringExpression || (expression instanceof ConstantExpression && expression.value instanceof String) )
            return thisMethodCall(expression, arguments)

        // e.g. 1 { }, 1.1 { }, (1 / 2) { }, m() { }, { -> ... } { }
        return callMethodCall(expression, arguments)
    }

    private Expression pathArgumentsElement(Expression caller, ArgumentListContext ctx) {
        final arguments = argumentList(ctx)
        return methodCall(caller, arguments)
    }

    private Expression pathIndexElement(Expression expression, IndexPropertyArgsContext ctx) {
        final elements = expressionList(ctx.expressionList())

        Expression index
        if( elements.size() > 1 ) {
            // e.g. a[1, 2]
            index = listX(elements)
            index.setWrapped(true)
        }
        else if( elements.first() instanceof SpreadExpression ) {
            // e.g. a[*[1, 2]]
            index = listX(elements)
            index.setWrapped(false)
        }
        else {
            // e.g. a[1]
            index = elements.first()
        }

        return indexX(expression, ast(index, ctx))
    }

    /// -- PRIMARY EXPRESSIONS

    private Expression primary(PrimaryContext ctx) {
        if( ctx instanceof IdentifierPrmrAltContext )
            return ast( varX(identifier(ctx.identifier())), ctx )

        if( ctx instanceof LiteralPrmrAltContext )
            return ast( literal(ctx.literal()), ctx )

        if( ctx instanceof GstringPrmrAltContext )
            return ast( gstring(ctx.gstring()), ctx )

        if( ctx instanceof NewPrmrAltContext )
            return ast( creator(ctx.creator()), ctx )

        if( ctx instanceof ParenPrmrAltContext )
            return ast( parExpression(ctx.parExpression()), ctx )

        if( ctx instanceof ClosurePrmrAltContext )
            return ast( closure(ctx.closure()), ctx )

        if( ctx instanceof ListPrmrAltContext )
            return ast( list(ctx.list()), ctx )

        if( ctx instanceof MapPrmrAltContext )
            return ast( map(ctx.map()), ctx )

        if( ctx instanceof BuiltInTypePrmrAltContext )
            return ast( builtInType(ctx.builtInType()), ctx )

        throw createParsingFailedException("Invalid Groovy expression: ${ctx.text}", ctx)
    }

    private Expression builtInType(BuiltInTypeContext ctx) {
        varX(ctx.text)
    }

    private String identifier(IdentifierContext ctx) {
        ctx.text
    }

    private String keywords(KeywordsContext ctx) {
        ctx.text
    }

    private Expression literal(LiteralContext ctx) {
        if( ctx instanceof IntegerLiteralAltContext )
            return integerLiteral( ctx )

        if( ctx instanceof FloatingPointLiteralAltContext )
            return floatingPointLiteral( ctx )

        if( ctx instanceof StringLiteralAltContext )
            return constX( stringLiteral(ctx.stringLiteral()) )

        if( ctx instanceof BooleanLiteralAltContext )
            return constX( ctx.text=='true' )

        if( ctx instanceof NullLiteralAltContext )
            return constX( null )

        throw createParsingFailedException("Invalid Groovy expression: ${ctx.text}", ctx)
    }

    private Expression integerLiteral(IntegerLiteralAltContext ctx) {
        Number num
        try {
            num = Numbers.parseInteger(ctx.text)
        } catch (Exception e) {
            numberFormatError = new Tuple2(ctx, e)
        }

        constX(num, true)
    }

    private Expression floatingPointLiteral(FloatingPointLiteralAltContext ctx) {
        Number num
        try {
            num = Numbers.parseDecimal(ctx.text)
        } catch (Exception e) {
            numberFormatError = new Tuple2(ctx, e)
        }

        constX(num, true)
    }

    private String stringLiteral(StringLiteralContext ctx) {
        stringLiteral(ctx.text)
    }

    private String stringLiteral(String text) {
        final startsWithSlash = text.startsWith(SLASH_STR)

        if( text.startsWith(TSQ_STR) || text.startsWith(TDQ_STR) ) {
            text = StringUtils.removeCR(text)
            text = StringUtils.trimQuotations(text, 3)
        }
        else if( text.startsWith(SQ_STR) || text.startsWith(DQ_STR) || startsWithSlash ) {
            // the slashy string can span rows, so we have to remove CR for it
            if( startsWithSlash )
                text = StringUtils.removeCR(text)
            text = StringUtils.trimQuotations(text, 1)
        }

        final slashyType = startsWithSlash
            ? StringUtils.SLASHY
            : StringUtils.NONE_SLASHY

        return StringUtils.replaceEscapes(text, slashyType)
    }

    private Expression gstring(GstringContext ctx) {
        final quote = gstringQuote(ctx.GStringBegin().text)
        final size = ctx.GStringPart().size() + 2

        final List<ConstantExpression> strings = new ArrayList(size)
        strings << gstringBegin(ctx.GStringBegin(), quote)
        for( final part : ctx.GStringPart() )
            strings << gstringPart(part, quote)
        strings << gstringEnd(ctx.GStringEnd(), quote)

        final verbatimText = stringLiteral(ctx.text)
        final values = ctx.gstringValue().collect( this.&gstringValue )
        new GStringExpression(verbatimText, strings, values)
    }

    private String gstringQuote(String text) {
        if( text.startsWith(TDQ_STR) )
            return TDQ_STR

        if( text.startsWith(DQ_STR) )
            return DQ_STR

        if( text.startsWith(SLASH_STR) )
            return SLASH_STR

        return String.valueOf(text.charAt(0))
    }

    private ConstantExpression gstringBegin(TerminalNode e, String quote) {
        final text = new StringBuilder(e.text)
        text.deleteCharAt(text.length() - 1)  // remove the trailing $
        text.append(quote)
        return ast( constX(stringLiteral(text.toString())), e )
    }

    private ConstantExpression gstringPart(TerminalNode e, String quote) {
        final text = new StringBuilder(e.text)
        text.deleteCharAt(text.length() - 1)  // remove the trailing $
        text.insert(0, quote).append(quote)
        return ast( constX(stringLiteral(text.toString())), e )
    }

    private ConstantExpression gstringEnd(TerminalNode e, String quote) {
        final text = new StringBuilder(e.text)
        text.insert(0, quote)
        return ast( constX(stringLiteral(text.toString())), e )
    }

    private Expression gstringValue(GstringValueContext ctx) {
        expression(ctx.expression())
    }

    private Expression creator(CreatorContext ctx) {
        final type = type(ctx.createdName())
        final arguments = argumentList(ctx.arguments().argumentList())
        ctorX(type, arguments)
    }

    private Expression parExpression(ParExpressionContext ctx) {
        final expression = expression(ctx.expression())
        expression.getNodeMetaData(INSIDE_PARENTHESES_LEVEL, k -> new AtomicInteger()).getAndAdd(1)
        return expression
    }

    private Expression closure(ClosureContext ctx) {
        final params = parameters(ctx.formalParameterList())
        final code = closureStatements(ctx.closureStatements())
        closureX(params, code)
    }

    private BlockStatement closureStatements(ClosureStatementsContext ctx) {
        final List<Statement> code = ctx
            ? ctx.statement().collect( this.&statement )
            : List<Statement>.of()
        ast( new BlockStatement(code, new VariableScope()), ctx )
    }

    private Expression list(ListContext ctx) {
        if( ctx.COMMA() && !ctx.expressionList() )
            throw createParsingFailedException("Empty list literal should not contain any comma(,)", ctx.COMMA())

        listX(expressionList(ctx.expressionList()))
    }

    private List<Expression> expressionList(ExpressionListContext ctx) {
        if( !ctx )
            return Collections.emptyList()
        
        ctx.expressionListElement().collect( this.&listElement )
    }

    private Expression listElement(ExpressionListElementContext ctx) {
        final element = expression(ctx.expression())
        ctx.MUL()
            ? ast( new SpreadExpression(element), ctx )
            : ast( element, ctx )
    }

    private Expression map(MapContext ctx) {
        if( !ctx.mapEntryList() )
            return new MapExpression()

        final entries = ctx.mapEntryList().mapEntry().collect( this.&mapEntry )
        mapX(entries)
    }

    private MapEntryExpression mapEntry(MapEntryContext ctx) {
        final value = expression(ctx.expression())
        final key = ctx.MUL()
            ? ast( new SpreadMapExpression(value), ctx )
            : mapEntryLabel(ctx.mapEntryLabel())
        ast( entryX(key, value), ctx )
    }

    private Expression mapEntryLabel(MapEntryLabelContext ctx) {
        if( ctx.keywords() )
            return ast( constX(keywords(ctx.keywords())), ctx )

        if( ctx.primary() ) {
            final expression = primary(ctx.primary())
            return expression instanceof VariableExpression && !isInsideParentheses(expression)
                ? ast( constX(((VariableExpression)expression).name), expression )
                : ast( expression, ctx )
        }

        throw createParsingFailedException("Unsupported map entry label: ${ctx.text}", ctx)
    }

    private Expression methodCall(Expression caller, Expression arguments) {
        // e.g. (obj.x)(), (obj.@x)()
        if( isInsideParentheses(caller) )
            return callMethodCall(caller, arguments)

        // e.g. obj.a(1, 2)
        if( caller instanceof PropertyExpression )
            return propMethodCall(caller, arguments)

        // e.g. m(), "$m"(), "m"()
        if( caller instanceof VariableExpression || caller instanceof GStringExpression || (caller instanceof ConstantExpression && caller.value instanceof String) )
            return thisMethodCall(caller, arguments)

        // e.g. 1(), 1.1(), ((int) 1 / 2)(1, 2), {a, b -> a + b }(1, 2), m()()
        return callMethodCall(caller, arguments)
    }

    private Expression propMethodCall(PropertyExpression caller, Expression arguments) {
        final result = callX( caller.objectExpression, caller.property, arguments )
        result.setImplicitThis(false)
        result.setSafe(caller.isSafe())
        result.setSpreadSafe(caller.isSpreadSafe())

        // method call obj*.m() -> safe=false and spreadSafe=true
        // property access obj*.p -> safe=true and spreadSafe=true
        if( caller.isSpreadSafe() )
            result.setSafe(false)

        return result
    }

    private Expression thisMethodCall(Expression caller, Expression arguments) {
        final object = varX('this')
        object.setColumnNumber(caller.getColumnNumber())
        object.setLineNumber(caller.getLineNumber())

        final name = caller instanceof VariableExpression
            ? constX(caller.text)
            : caller

        return callX( object, name, arguments )
    }

    private Expression callMethodCall(Expression caller, Expression arguments) {
        final call = callX(caller, CALL_STR, arguments)
        call.setImplicitThis(false)
        return call
    }

    private Expression argumentList(ArgumentListContext ctx) {
        if( !ctx )
            return new ArgumentListExpression()

        final List<Expression> arguments = []
        final List<MapEntryExpression> opts = []

        for( final ctx1 : ctx.argumentListElement() ) {
            if( ctx1.expressionListElement() )
                arguments << listElement(ctx1.expressionListElement())

            else if( ctx1.namedArg() )
                opts << namedArg(ctx1.namedArg())

            else
                throw createParsingFailedException("Invalid Groovy method argument: ${ctx.text}", ctx)
        }

        // TODO: validate duplicate named arguments ?
        // TODO: only named arguments -> TupleExpression ?
        if( opts )
            arguments.push( ast(mapX(opts), ctx) )

        return ast( args(arguments), ctx )
    }

    private MapEntryExpression namedArg(NamedArgContext ctx) {
        final value = expression(ctx.expression())
        final key = ctx.MUL()
            ? new SpreadMapExpression(value)
            : namedArgLabel(ctx.namedArgLabel())
        new MapEntryExpression(key, value)
    }

    private Expression namedArgLabel(NamedArgLabelContext ctx) {
        if( ctx.keywords() )
            return constX(keywords(ctx.keywords()))

        if( ctx.identifier() )
            return constX(identifier(ctx.identifier()))

        if( ctx.literal() )
            return literal(ctx.literal())

        if( ctx.gstring() )
            return gstring(ctx.gstring())

        throw createParsingFailedException("Invalid Groovy method named argument: ${ctx.text}", ctx)
    }

    /// MISCELLANEOUS

    private Parameter[] parameters(FormalParameterListContext ctx) {
        // NOTE: implicit `it` is not allowed
        if( !ctx )
            return null

        for( int i = 0, n = ctx.formalParameter().size(); i < n - 1; i += 1 ) {
            final ctx1 = ctx.formalParameter(i)
            if( ctx1.ELLIPSIS() )
                throw createParsingFailedException("The var-arg parameter must be the last parameter", ctx1)
        }

        final params = ctx.formalParameter().collect( this.&parameter )
        for( int n = params.size(), i = n - 1; i >= 0; i -= 1 ) {
            final param = params[i]
            for( final other : params ) {
                if( other == param )
                    continue
                if( other.name == param.name )
                    throw createParsingFailedException("Duplicated parameter '${param.name}' found", param)
            }
        }

        return params as Parameter[]
    }

    private Parameter parameter(FormalParameterContext ctx) {
        ClassNode type = type(ctx.type())
        if( ctx.ELLIPSIS() ) {
            type = type.makeArray()
            type = ctx.type()
                ? ast( type, ctx.type(), ast(constX('...'), ctx.ELLIPSIS()) )
                : ast( type, ctx.ELLIPSIS() )
        }

        final name = identifier(ctx.identifier())
        final defaultValue = ctx.expression()
            ? expression(ctx.expression())
            : null
        ast( param(type, name, defaultValue), ctx )
    }

    private Token token(ParserToken token, int cardinality = 1) {
        final text = cardinality == 1
            ? token.text
            : token.text * cardinality
        final type = token.type == RANGE_EXCLUSIVE_RIGHT || token.type == RANGE_INCLUSIVE
            ? Types.RANGE_OPERATOR
            : Types.lookup(text, Types.ANY)
        new Token( type, text, token.getLine(), token.getCharPositionInLine() + 1 )
    }

    private ClassNode type(CreatedNameContext ctx) {
        if( ctx.qualifiedClassName() ) {
            final classNode = type(ctx.qualifiedClassName())
            if( ctx.typeArgumentsOrDiamond() )
                classNode.setGenericsTypes( typeArguments(ctx.typeArgumentsOrDiamond()) )
            return ast( classNode, ctx )
        }

        if( ctx.primitiveType() )
            return ast( type(ctx.primitiveType()), ctx )

        throw createParsingFailedException("Unrecognized created name: ${ctx.text}", ctx)
    }

    private ClassNode type(PrimitiveTypeContext ctx) {
        ClassHelper.make(ctx.text).getPlainNodeReference(false)
    }

    private ClassNode type(QualifiedClassNameContext ctx, boolean allowProxy=true) {
        final classNode = ClassHelper.make(ctx.text)

        if( classNode.isUsingGenerics() && allowProxy ) {
            final proxy = ClassHelper.makeWithoutCaching(classNode.name)
            proxy.setRedirect(classNode)
            return proxy
        }

        return classNode
    }

    private ClassNode type(TypeContext ctx, boolean allowProxy=true) {
        if( !ctx )
            return ClassHelper.dynamicType()

        if( ctx.qualifiedClassName() ) {
            final classNode = type(ctx.qualifiedClassName(), allowProxy)
            if( ctx.typeArguments() )
                classNode.setGenericsTypes( typeArguments(ctx.typeArguments()) )
            return ast( classNode, ctx )
        }

        if( ctx.primitiveType() )
            return ast( type(ctx.primitiveType()), ctx )

        throw createParsingFailedException("Unrecognized type: ${ctx.text}", ctx)
    }

    private GenericsType[] typeArguments(TypeArgumentsOrDiamondContext ctx) {
        ctx.typeArguments()
            ? typeArguments(ctx.typeArguments())
            : GenericsType.EMPTY_ARRAY
    }

    private GenericsType[] typeArguments(TypeArgumentsContext ctx) {
        ctx.type().collect( this.&genericsType ) as GenericsType[]
    }

    private GenericsType genericsType(TypeContext ctx) {
        ast( new GenericsType(type(ctx)), ctx )
    }

    /// HELPERS

    private boolean isInsideParentheses(NodeMetaDataHandler nodeMetaDataHandler) {
        Number insideParenLevel = nodeMetaDataHandler.getNodeMetaData(INSIDE_PARENTHESES_LEVEL)
        return insideParenLevel != null && insideParenLevel.intValue() > 0
    }

    private CompilationFailedException createParsingFailedException(String msg, ParserRuleContext ctx) {
        return createParsingFailedException(
            new SyntaxException(msg,
                ctx.start.getLine(),
                ctx.start.getCharPositionInLine() + 1,
                ctx.stop.getLine(),
                ctx.stop.getCharPositionInLine() + 1 + ctx.stop.getText().length()))
    }

    private CompilationFailedException createParsingFailedException(String msg, Tuple2<Integer,Integer> start, Tuple2<Integer,Integer> end) {
        return createParsingFailedException(
            new SyntaxException(msg,
                start.getV1(),
                start.getV2(),
                end.getV1(),
                end.getV2()))
    }

    private CompilationFailedException createParsingFailedException(String msg, ASTNode node) {
        return createParsingFailedException(
            new SyntaxException(msg,
                node.getLineNumber(),
                node.getColumnNumber(),
                node.getLastLineNumber(),
                node.getLastColumnNumber()))
    }

    private CompilationFailedException createParsingFailedException(String msg, TerminalNode node) {
        return createParsingFailedException(msg, node.getSymbol())
    }

    private CompilationFailedException createParsingFailedException(String msg, ParserToken token) {
        return createParsingFailedException(
            new SyntaxException(msg,
                token.getLine(),
                token.getCharPositionInLine() + 1,
                token.getLine(),
                token.getCharPositionInLine() + 1 + token.getText().length()))
    }

    private CompilationFailedException createParsingFailedException(Throwable t) {
        if( t instanceof SyntaxException )
            this.collectSyntaxError(t)

        else if( t instanceof GroovySyntaxError )
            this.collectSyntaxError(
                    new SyntaxException(
                            t.getMessage(),
                            t,
                            t.getLine(),
                            t.getColumn()))

        else if( t instanceof Exception )
            this.collectException(t)

        return new CompilationFailedException(
                CompilePhase.PARSING.getPhaseNumber(),
                this.sourceUnit,
                t)
    }

    private void collectSyntaxError(SyntaxException e) {
        sourceUnit.getErrorCollector().addFatalError(new SyntaxErrorMessage(e, sourceUnit))
    }

    private void collectException(Exception e) {
        sourceUnit.getErrorCollector().addException(e, this.sourceUnit)
    }

    private void removeErrorListeners() {
        lexer.removeErrorListeners()
        parser.removeErrorListeners()
    }

    private void addErrorListeners() {
        lexer.removeErrorListeners()
        lexer.addErrorListener(createANTLRErrorListener())

        parser.removeErrorListeners()
        parser.addErrorListener(createANTLRErrorListener())
    }

    private ANTLRErrorListener createANTLRErrorListener() {
        return new ANTLRErrorListener() {
            @Override
            void syntaxError(Recognizer recognizer, Object offendingSymbol, int line, int charPositionInLine, String msg, RecognitionException e) {
                collectSyntaxError(new SyntaxException(msg, line, charPositionInLine + 1))
            }

            @Override
            void reportAmbiguity(Parser recognizer, DFA dfa, int startIndex, int stopIndex, boolean exact, BitSet ambigAlts, ATNConfigSet configs) {}

            @Override
            void reportAttemptingFullContext(Parser recognizer, DFA dfa, int startIndex, int stopIndex, BitSet conflictingAlts, ATNConfigSet configs) {}

            @Override
            void reportContextSensitivity(Parser recognizer, DFA dfa, int startIndex, int stopIndex, int prediction, ATNConfigSet configs) {}
        }
    }

    private static final String CALL_STR = 'call'
    private static final String SLASH_STR = '/'
    private static final String TDQ_STR = '"""'
    private static final String TSQ_STR = "'''"
    private static final String SQ_STR = "'"
    private static final String DQ_STR = '"'

    private static final String INSIDE_PARENTHESES_LEVEL = "_INSIDE_PARENTHESES_LEVEL"

}
