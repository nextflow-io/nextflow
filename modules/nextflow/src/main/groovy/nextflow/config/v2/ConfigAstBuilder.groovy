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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.antlr.ConfigLexer
import nextflow.antlr.ConfigParser
import nextflow.antlr.DescriptiveErrorStrategy
import org.antlr.v4.runtime.CharStream
import org.antlr.v4.runtime.CharStreams
import org.antlr.v4.runtime.CommonTokenStream
import org.antlr.v4.runtime.ParserRuleContext
import org.antlr.v4.runtime.Token as ParserToken
import org.antlr.v4.runtime.atn.PredictionMode
import org.antlr.v4.runtime.misc.ParseCancellationException
import org.antlr.v4.runtime.tree.TerminalNode
import org.apache.groovy.parser.antlr4.GroovySyntaxError
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassHelper
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.ModuleNode
import org.codehaus.groovy.ast.Parameter
import org.codehaus.groovy.ast.VariableScope
import org.codehaus.groovy.ast.expr.ArgumentListExpression
import org.codehaus.groovy.ast.expr.AttributeExpression
import org.codehaus.groovy.ast.expr.BinaryExpression
import org.codehaus.groovy.ast.expr.BitwiseNegationExpression
import org.codehaus.groovy.ast.expr.BooleanExpression
import org.codehaus.groovy.ast.expr.CastExpression
import org.codehaus.groovy.ast.expr.ClassExpression
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.ConstantExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.GStringExpression
import org.codehaus.groovy.ast.expr.LambdaExpression
import org.codehaus.groovy.ast.expr.ListExpression
import org.codehaus.groovy.ast.expr.MapExpression
import org.codehaus.groovy.ast.expr.MapEntryExpression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.NotExpression
import org.codehaus.groovy.ast.expr.PostfixExpression
import org.codehaus.groovy.ast.expr.PrefixExpression
import org.codehaus.groovy.ast.expr.PropertyExpression
import org.codehaus.groovy.ast.expr.SpreadExpression
import org.codehaus.groovy.ast.expr.SpreadMapExpression
import org.codehaus.groovy.ast.expr.TernaryExpression
import org.codehaus.groovy.ast.expr.TupleExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.ast.expr.UnaryMinusExpression
import org.codehaus.groovy.ast.expr.UnaryPlusExpression
import org.codehaus.groovy.ast.stmt.AssertStatement
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.EmptyStatement
import org.codehaus.groovy.ast.stmt.ReturnStatement
import org.codehaus.groovy.ast.stmt.Statement
import org.codehaus.groovy.ast.tools.GeneralUtils
import org.codehaus.groovy.control.CompilationFailedException
import org.codehaus.groovy.control.CompilePhase
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.control.messages.SyntaxErrorMessage
import org.codehaus.groovy.syntax.Numbers
import org.codehaus.groovy.syntax.SyntaxException
import org.codehaus.groovy.syntax.Token
import org.codehaus.groovy.syntax.Types

import static nextflow.antlr.ConfigParser.*
import static org.codehaus.groovy.ast.tools.GeneralUtils.args as argsX
import static org.codehaus.groovy.ast.tools.GeneralUtils.callX
import static org.codehaus.groovy.ast.tools.GeneralUtils.callThisX
import static org.codehaus.groovy.ast.tools.GeneralUtils.closureX
import static org.codehaus.groovy.ast.tools.GeneralUtils.constX
import static org.codehaus.groovy.ast.tools.GeneralUtils.ctorX
import static org.codehaus.groovy.ast.tools.GeneralUtils.declX
import static org.codehaus.groovy.ast.tools.GeneralUtils.stmt
import static org.codehaus.groovy.ast.tools.GeneralUtils.varX

@Slf4j
@CompileStatic
class ConfigAstBuilder {

    private SourceUnit sourceUnit
    private ModuleNode moduleNode
    private ConfigLexer lexer
    private ConfigParser parser

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
        // TODO: configure script class node ?
        // TODO: check number format error
        return moduleNode
    }

    private Statement configStatement(ConfigStatementContext ctx) {
        if( ctx instanceof ConfigIncludeStmtAltContext )
            return configInclude(ctx.configInclude())

        if( ctx instanceof ConfigAssignmentStmtAltContext )
            return configAssignment(ctx.configAssignment())

        if( ctx instanceof ConfigBlockStmtAltContext )
            return configBlock(ctx.configBlock())

        throw createParsingFailedException("Invalid config statement: ${ctx.text}", ctx)
    }

    private Statement configInclude(ConfigIncludeContext ctx) {
        final source = expression(ctx.expression())
        final include = callThisX('includeConfig', argsX(source))
        stmt(include)
    }

    private Statement configAssignment(ConfigAssignmentContext ctx) {
        final names = new ListExpression( ctx.configPathExpression().identifier().collect( ctx1 -> (Expression)constX(identifier(ctx1)) ) )
        final right = expression(ctx.expression())
        stmt(callThisX('assign', argsX([names, right])))
    }

    private Statement configBlock(ConfigBlockContext ctx) {
        final name = ctx.identifier()
            ? constX(identifier(ctx.identifier()))
            : constX(stringLiteral(ctx.stringLiteral()))
        final statements = ctx.configBlockStatement().collect( ctx1 -> configBlockStatement(ctx1) )
        final closure = closureX(new BlockStatement(statements, new VariableScope()))
        stmt(callThisX('block', argsX([name, closure])))
    }

    private Statement configBlockStatement(ConfigBlockStatementContext ctx) {
        if( ctx instanceof ConfigIncludeBlockStmtAltContext )
            return configInclude(ctx.configInclude())

        if( ctx instanceof ConfigAssignmentBlockStmtAltContext )
            return configAssignment(ctx.configAssignment())

        if( ctx instanceof ConfigBlockBlockStmtAltContext )
            return configBlock(ctx.configBlock())

        if( ctx instanceof ConfigSelectorBlockStmtAltContext )
            return configSelector(ctx.configSelector())

        throw createParsingFailedException("Invalid statement in config block: ${ctx.text}", ctx)
    }

    private Statement configSelector(ConfigSelectorContext ctx) {
        final kind = ctx.kind.text
        final target = configSelectorTarget(ctx.target)
        final statements = ctx.configAssignment().collect( ctx1 -> configAssignment(ctx1) )
        final closure = closureX(new BlockStatement(statements, new VariableScope()))
        stmt(callThisX(kind, argsX([target, closure])))
    }

    private Expression configSelectorTarget(ConfigSelectorTargetContext ctx) {
        ctx.identifier()
            ? constX(identifier(ctx.identifier()))
            : constX(stringLiteral(ctx.stringLiteral()))
    }

    /// GROOVY STATEMENTS

    private List<Statement> statements(List<StatementContext> stmts) {
        stmts.collect( ctx -> statement(ctx) )
    }

    private Statement statement(StatementContext ctx) {
        if( ctx instanceof BlockStmtAltContext )
            return block(ctx.block())

        if( ctx instanceof ReturnStmtAltContext )
            return returnStatement(ctx.expression())

        if( ctx instanceof AssertStmtAltContext )
            return assertStatement(ctx.assertStatement())

        if( ctx instanceof VariableDeclarationStmtAltContext )
            return variableDeclaration(ctx.variableDeclaration())

        if( ctx instanceof ExpressionStmtAltContext )
            return expressionStatement(ctx.expressionStatement())

        if( ctx instanceof EmptyStmtAltContext )
            return new EmptyStatement()

        throw createParsingFailedException("Invalid Groovy statement: ${ctx.text}", ctx)
    }

    private BlockStatement block(BlockContext ctx) {
        block(ctx.blockStatements())
    }

    private BlockStatement block(BlockStatementsContext ctx) {
        final code = ctx
            ? statements(ctx.statement())
            : [ new EmptyStatement() as Statement ]
        new BlockStatement(code, new VariableScope())
    }

    private ReturnStatement returnStatement(ExpressionContext ctx) {
        final result = ctx
            ? expression(ctx)
            : ConstantExpression.EMPTY_EXPRESSION
        new ReturnStatement(result)
    }

    private AssertStatement assertStatement(AssertStatementContext ctx) {
        final condition = new BooleanExpression(expression(ctx.condition))
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
                varX(name, type)
            }
            final target = variables.size() > 1
                ? new ArgumentListExpression(variables as List<Expression>)
                : variables.first()
            final initializer = expression(ctx.initializer)
            return stmt(declX(target, initializer))
        }
        else {
            // single assignment
            final type = type(ctx.type())
            final decl = ctx.variableDeclarator()
            final name = identifier(decl.identifier())
            final target = varX(name, type)
            final initializer = expression(decl.initializer)
            return stmt(declX(target, initializer))
        }
    }

    private Statement expressionStatement(ExpressionStatementContext ctx) {
        ctx.argumentList()
            ? stmt(methodCall(ctx.expression(), ctx.argumentList()))
            : stmt(expression(ctx.expression()))
    }

    /// GROOVY EXPRESSIONS

    private Expression expression(ExpressionContext ctx) {
        if( ctx instanceof AddExprAltContext )
            return binary(ctx.left, ctx.op, ctx.right)

        if( ctx instanceof AndExprAltContext )
            return binary(ctx.left, ctx.op, ctx.right)

        if( ctx instanceof AssignmentExprAltContext )
            return binary(ctx.left, ctx.op, ctx.right)

        if( ctx instanceof CastExprAltContext ) {
            final type = type(ctx.castParExpression().type())
            final operand = castOperand(ctx.castOperandExpression())
            return new CastExpression(type, operand)
        }

        if( ctx instanceof ConditionalExprAltContext )
            return ternary(ctx.condition, ctx.tb, ctx.fb)

        if( ctx instanceof EqualityExprAltContext )
            return binary(ctx.left, ctx.op, ctx.right)

        if( ctx instanceof ExclusiveOrExprAltContext )
            return binary(ctx.left, ctx.op, ctx.right)

        if( ctx instanceof InclusiveOrExprAltContext )
            return binary(ctx.left, ctx.op, ctx.right)

        if( ctx instanceof LogicalAndExprAltContext )
            return binary(ctx.left, ctx.op, ctx.right)

        if( ctx instanceof LogicalOrExprAltContext )
            return binary(ctx.left, ctx.op, ctx.right)

        if( ctx instanceof MultDivExprAltContext )
            return binary(ctx.left, ctx.op, ctx.right)

        if( ctx instanceof MultipleAssignmentExprAltContext ) {
            final vars = ctx.variableNames().identifier().collect( ctx1 -> varX(identifier(ctx1)) )
            final left = new TupleExpression(vars as List<Expression>)
            final right = expression(ctx.right)
            return new BinaryExpression(left, token(ctx.op), right)
        }

        if( ctx instanceof PathExprAltContext )
            return path(ctx.pathExpression())

        if( ctx instanceof PostfixExprAltContext )
            return postfix(ctx.pathExpression(), ctx.op)

        if( ctx instanceof PowerExprAltContext )
            return binary(ctx.left, ctx.op, ctx.right)

        if( ctx instanceof PrefixExprAltContext )
            return prefix(expression(ctx.expression()), ctx.op)

        if( ctx instanceof RegexExprAltContext )
            return binary(ctx.left, ctx.op, ctx.right)

        if( ctx instanceof RelationalExprAltContext ) {
            final right = ctx.type()
                ? new ClassExpression(type(ctx.type()))
                : expression(ctx.right)
            return binary(ctx.left, ctx.op, right)
        }

        if( ctx instanceof ShiftExprAltContext ) {
            final op = ctx.dlOp ?: ctx.tgOp ?: ctx.dgOp ?: ctx.riOp ?: ctx.reOp
            return binary(ctx.left, op, ctx.right)
        }

        if( ctx instanceof UnaryAddExprAltContext )
            return unaryAdd(expression(ctx.expression()), ctx.op, ctx)

        if( ctx instanceof UnaryNotExprAltContext )
            return unaryNot(expression(ctx.expression()), ctx.op, ctx)

        throw createParsingFailedException("Invalid Groovy expression: ${ctx.text}", ctx)
    }

    private BinaryExpression binary(ExpressionContext left, ParserToken op, ExpressionContext right) {
        new BinaryExpression(expression(left), token(op), expression(right))
    }

    private BinaryExpression binary(ExpressionContext left, ParserToken op, Expression right) {
        new BinaryExpression(expression(left), token(op), right)
    }

    private Expression castOperand(CastOperandExpressionContext ctx) {
        if( ctx instanceof CastCastExprAltContext ) {
            final type = type(ctx.castParExpression().type())
            final operand = castOperand(ctx.castOperandExpression())
            return new CastExpression(type, operand)
        }

        if( ctx instanceof PathCastExprAltContext )
            return path(ctx.pathExpression())

        if( ctx instanceof PostfixCastExprAltContext )
            return postfix(ctx.pathExpression(), ctx.op)

        if( ctx instanceof PrefixCastExprAltContext )
            return prefix(castOperand(ctx.castOperandExpression()), ctx.op)

        if( ctx instanceof UnaryAddCastExprAltContext )
            return unaryAdd(castOperand(ctx.castOperandExpression()), ctx.op, ctx)

        if( ctx instanceof UnaryNotCastExprAltContext )
            return unaryNot(castOperand(ctx.castOperandExpression()), ctx.op, ctx)

        throw createParsingFailedException("Invalid Groovy expression: ${ctx.text}", ctx)
    }

    private PostfixExpression postfix(PathExpressionContext ctx, ParserToken op) {
        new PostfixExpression(path(ctx), token(op))
    }

    private PrefixExpression prefix(Expression expression, ParserToken op) {
        new PrefixExpression(token(op), expression)
    }

    private TernaryExpression ternary(ExpressionContext cond, ExpressionContext tb, ExpressionContext fb) {
        final condition = new BooleanExpression(expression(cond))
        final trueExpression = tb ? expression(tb) : null
        final falseExpression = expression(fb)
        new TernaryExpression(condition, trueExpression, falseExpression)
    }

    private Expression unaryAdd(Expression expression, ParserToken op, ParserRuleContext ctx) {
        if( op.type == ConfigParser.ADD )
            return isNonStringConstantOutsideParentheses(expression)
                ? expression
                : new UnaryPlusExpression(expression)

        if( op.type == ConfigParser.SUB )
            return isNonStringConstantOutsideParentheses(expression)
                ? constX(((ConstantExpression)expression).value, true)
                : new UnaryMinusExpression(expression)

        throw createParsingFailedException("Unsupported unary expression: ${ctx.text}", ctx)
    }

    private boolean isNonStringConstantOutsideParentheses(Expression expression) {
        return expression instanceof ConstantExpression
                && expression.value !instanceof String
                // && !isInsideParentheses(expression)
    }

    // private boolean isInsideParentheses(NodeMetaDataHandler nodeMetaDataHandler) {
    //     final insideParenLevel = nodeMetaDataHandler.getNodeMetaData(INSIDE_PARENTHESES_LEVEL)
    //     return insideParenLevel != null && insideParenLevel.intValue() > 0
    // }

    private Expression unaryNot(Expression expression, ParserToken op, ParserRuleContext ctx) {
        if( op.type == ConfigParser.NOT )
            return new NotExpression(expression)

        if( op.type == ConfigParser.BITNOT )
            return new BitwiseNegationExpression(expression)

        throw createParsingFailedException("Unsupported unary expression: ${ctx.text}", ctx)
    }

    /// -- PATH EXPRESSIONS

    private Expression path(PathExpressionContext ctx) {
        try {
            return ctx.pathElement().inject(primary(ctx.primary()), (acc, el) -> pathElement(acc, el))
        }
        catch( IllegalStateException e ) {
            throw createParsingFailedException("Invalid Groovy expression: ${ctx.text}", ctx)
        }
    }

    private Expression pathElement(Expression expression, PathElementContext ctx) {
        if( ctx instanceof PropertyPathExprAltContext ) {
            final prop =
                ctx.keywords() ? ctx.keywords().text :
                ctx.identifier() ? identifier(ctx.identifier()) :
                /* ctx.stringLiteral() */ stringLiteral(ctx.stringLiteral())
            final safe = ctx.SAFE_DOT() != null || ctx.SPREAD_DOT() != null
            final result = new PropertyExpression(expression, constX(prop), safe)
            if( ctx.SPREAD_DOT() )
                result.setSpreadSafe(true)
            return result
        }

        if( ctx instanceof ClosurePathExprAltContext ) {
            if( expression instanceof MethodCallExpression ) {
                // append closure to method call arguments
                final methodCall = (MethodCallExpression)expression
                final arguments = (ArgumentListExpression)methodCall.arguments
                arguments.addExpression( closure(ctx.closure()) )
                return methodCall
            }
            else {
                // create method call with single closure argument
                final parts = methodObject(expression)
                final closure = closure(ctx.closure())
                return callX(parts.first, parts.second, argsX(closure))
            }
        }

        if( ctx instanceof ArgumentsPathExprAltContext ) {
            final parts = methodObject(expression)
            final arguments = argumentList(ctx.arguments().argumentList())
            return callX(parts.first, parts.second, arguments)
        }

        if( ctx instanceof ListElementPathExprAltContext ) {
            final op = new Token(Types.LEFT_SQUARE_BRACKET, '[', -1, -1)
            final elements = expressionList(ctx.expressionList())
            final arg = elements.size() > 1
                ? new ListExpression(elements)
                : elements.first()
            final result = new BinaryExpression(expression, op, arg)
            if( ctx.QUESTION() )
                result.setSafe(true)
            return result
        }

        throw new IllegalStateException()
    }

    /// -- PRIMARY EXPRESSIONS

    private Expression primary(PrimaryContext ctx) {
        if( ctx instanceof IdentifierPrmrAltContext )
            return varX(identifier(ctx.identifier()))

        if( ctx instanceof LiteralPrmrAltContext )
            return constant(ctx.literal())

        if( ctx instanceof GstringPrmrAltContext )
            return gstring(ctx.gstring())

        if( ctx instanceof NewPrmrAltContext )
            return creator(ctx.creator())

        if( ctx instanceof ParenPrmrAltContext )
            return expression(ctx.parExpression().expression())

        if( ctx instanceof ClosurePrmrAltContext )
            return closure(ctx.closure())

        if( ctx instanceof ListPrmrAltContext )
            return list(ctx.list())

        if( ctx instanceof MapPrmrAltContext )
            return map(ctx.map())

        if( ctx instanceof BuiltInTypePrmrAltContext )
            return new ClassExpression(type(ctx.builtInType()))

        throw createParsingFailedException("Invalid Groovy expression: ${ctx.text}", ctx)
    }

    private String identifier(IdentifierContext ctx) {
        ctx.text
    }

    private ConstantExpression constant(LiteralContext ctx) {
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

    private ConstantExpression integerLiteral(IntegerLiteralAltContext ctx) {
        Number num
        try {
            num = Numbers.parseInteger(ctx.text)
        } catch (Exception e) {
            // this.numberFormatError = tuple(ctx, e)
        }

        constX(num, true)
    }

    private ConstantExpression floatingPointLiteral(FloatingPointLiteralAltContext ctx) {
        Number num
        try {
            num = Numbers.parseDecimal(ctx.text)
        } catch (Exception e) {
            // this.numberFormatError = tuple(ctx, e)
        }

        constX(num, true)
    }

    private String stringLiteral(StringLiteralContext ctx) {
        unquote(ctx.text)
    }

    private GStringExpression gstring(GstringContext ctx) {
        final strings = gstringParts(ctx)
        final values = ctx.gstringValue().collect( this.&gstringValue )
        new GStringExpression(unquote(ctx.text), strings, values)
    }

    private List<ConstantExpression> gstringParts(GstringContext ctx) {
        final strings = []

        final begin = ctx.GStringBegin().text
        strings.add(begin.startsWith('"""') ? begin[3..<-1] : begin[1..<-1])

        for( final part : ctx.GStringPart() )
            strings.add(part.text[0..<-1])

        final end = ctx.GStringEnd().text
        strings.add(begin.startsWith('"""') ? end[0..<-3] : end[0..<-1])

        return strings.collect( str -> constX(str) )
    }

    private Expression gstringValue(GstringValueContext ctx) {
        expression(ctx.expression())
    }

    private Expression creator(CreatorContext ctx) {
        final type = type(ctx.createdName())
        final arguments = argumentList(ctx.arguments().argumentList())
        ctorX(type, arguments)
    }

    private Expression closure(ClosureContext ctx) {
        final params = parameters(ctx.formalParameterList())
        final code = block(ctx.blockStatements())
        final closure = closureX(params, code)
        final source = constX(ctx.text)
        ctorX(type0('ClosureWithSource'), argsX([closure, source]))
        // source
    }

    private ListExpression list(ListContext ctx) {
        if( ctx.COMMA() && !ctx.expressionList() )
            throw createParsingFailedException("Empty list literal should not contain any comma(,)", ctx.COMMA())

        new ListExpression(expressionList(ctx.expressionList()))
    }

    private List<Expression> expressionList(ExpressionListContext ctx) {
        if( !ctx )
            return Collections.emptyList()
        
        ctx.expressionListElement().collect( this.&listElement )
    }

    private Expression listElement(ExpressionListElementContext ctx) {
        final element = expression(ctx.expression())
        ctx.MUL()
            ? new SpreadExpression(element)
            : element
    }

    private MapExpression map(MapContext ctx) {
        if( !ctx.mapEntryList() )
            return new MapExpression()

        final entries = ctx.mapEntryList().mapEntry().collect( this.&mapEntry )
        new MapExpression(entries)
    }

    private MapEntryExpression mapEntry(MapEntryContext ctx) {
        final value = expression(ctx.expression())
        final key = ctx.MUL()
            ? new SpreadMapExpression(value)
            : mapEntryLabel(ctx.mapEntryLabel())
        new MapEntryExpression(key, value)
    }

    private Expression mapEntryLabel(MapEntryLabelContext ctx) {
        if( ctx.keywords() )
            return constX(ctx.text)

        if( ctx.primary() )
            return primary(ctx.primary())

        throw createParsingFailedException("Unsupported map entry label: ${ctx.text}", ctx)
    }

    private MethodCallExpression methodCall(ExpressionContext method, ArgumentListContext args) {
        final parts = methodObject(expression(method))
        final arguments = argumentList(args)
        callX(parts.first, parts.second, arguments)
    }

    private Tuple2<Expression,Expression> methodObject(Expression expression) {
        if( expression instanceof PropertyExpression )
            return new Tuple2(expression.objectExpression, expression.property)

        if( expression instanceof VariableExpression )
            return new Tuple2(varX('this'), constX(expression.text))

        return new Tuple2(expression, constX('call'))
    }

    private Expression argumentList(ArgumentListContext ctx) {
        if( !ctx )
            return new ArgumentListExpression()

        final List<Expression> args = []
        final List<MapEntryExpression> opts = []

        for( final ctx1 : ctx.argumentListElement() ) {
            if( ctx1.expressionListElement() )
                args << listElement(ctx1.expressionListElement())

            else if( ctx1.namedArg() )
                opts << namedArg(ctx1.namedArg())

            else
                throw createParsingFailedException("Invalid Groovy method argument: ${ctx.text}", ctx)
        }

        if( opts )
            args.push(new MapExpression(opts))

        return new ArgumentListExpression(args)
    }

    private MapEntryExpression namedArg(NamedArgContext ctx) {
        final value = expression(ctx.expression())
        final key = ctx.MUL()
            ? new SpreadMapExpression(value)
            : namedArgLabel(ctx.namedArgLabel())
        new MapEntryExpression(key, value)
    }

    private Expression namedArgLabel(NamedArgLabelContext ctx) {
        if( ctx.keywords() || ctx.identifier() )
            return constX(ctx.text)

        if( ctx.literal() )
            return constant(ctx.literal())

        if( ctx.gstring() )
            return gstring(ctx.gstring())

        throw createParsingFailedException("Invalid Groovy method named argument: ${ctx.text}", ctx)
    }

    /// MISCELLANEOUS

    private Parameter[] parameters(FormalParameterListContext ctx) {
        final result = ctx
            ? ctx.formalParameter().collect( this.&parameter )
            : []
        return result as Parameter[]
    }

    private Parameter parameter(FormalParameterContext ctx) {
        final type = type(ctx.type())
        final name = identifier(ctx.identifier())
        final defaultValue = ctx.expression() ? expression(ctx.expression()) : null
        new Parameter(type, name, defaultValue)
    }

    private Token token(ParserToken t) {
        new Token(Types.lookupSymbol(t.text), t.text, -1, -1)
    }

    private ClassNode type(ParserRuleContext ctx) {
        ctx
            ? type0(ctx.text)
            : ClassHelper.DYNAMIC_TYPE
    }

    private ClassNode type0(String text) {
        ClassHelper.makeWithoutCaching(text)
    }

    /// HELPERS

    private String unquote(String text) {
        if( !text )
            return null
        else if( text.startsWith('"""') || text.startsWith("'''") )
            return text[3..<-3]
        else
            return text[1..<-1]
    }

    private CompilationFailedException createParsingFailedException(String msg, ParserRuleContext ctx) {
        return createParsingFailedException(
            new SyntaxException(msg,
                ctx.start.getLine(),
                ctx.start.getCharPositionInLine() + 1,
                ctx.stop.getLine(),
                ctx.stop.getCharPositionInLine() + 1 + ctx.stop.getText().length()))
    }

    private CompilationFailedException createParsingFailedException(String msg, Tuple2<Integer, Integer> start, Tuple2<Integer, Integer> end) {
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
        // TODO: missing ANTLRErrorListener::reportAmbiguity()

        // lexer.removeErrorListeners()
        // lexer.addErrorListener(this.createANTLRErrorListener())

        // parser.removeErrorListeners()
        // parser.addErrorListener(this.createANTLRErrorListener())
    }

}
