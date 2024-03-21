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

package nextflow.config

import groovy.transform.CompileStatic
import nextflow.antlr.ConfigLexer
import nextflow.antlr.ConfigParser
import org.antlr.v4.runtime.BailErrorStrategy
import org.antlr.v4.runtime.CharStream
import org.antlr.v4.runtime.CharStreams
import org.antlr.v4.runtime.CommonTokenStream
import org.antlr.v4.runtime.ParserRuleContext
import org.antlr.v4.runtime.Token as ParserToken
import org.antlr.v4.runtime.tree.Tree
import org.codehaus.groovy.ast.ClassHelper
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.Parameter
import org.codehaus.groovy.ast.VariableScope
import org.codehaus.groovy.ast.expr.ArgumentListExpression
import org.codehaus.groovy.ast.expr.AttributeExpression
import org.codehaus.groovy.ast.expr.BinaryExpression
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
import org.codehaus.groovy.ast.expr.PostfixExpression
import org.codehaus.groovy.ast.expr.PrefixExpression
import org.codehaus.groovy.ast.expr.PropertyExpression
import org.codehaus.groovy.ast.expr.SpreadExpression
import org.codehaus.groovy.ast.expr.SpreadMapExpression
import org.codehaus.groovy.ast.expr.TernaryExpression
import org.codehaus.groovy.ast.expr.TupleExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.ast.stmt.AssertStatement
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.EmptyStatement
import org.codehaus.groovy.ast.stmt.ReturnStatement
import org.codehaus.groovy.ast.stmt.Statement
import org.codehaus.groovy.ast.tools.GeneralUtils
import org.codehaus.groovy.syntax.Token
import org.codehaus.groovy.syntax.Types

import static nextflow.antlr.ConfigParser.*
import static org.codehaus.groovy.ast.tools.GeneralUtils.args as argsX
import static org.codehaus.groovy.ast.tools.GeneralUtils.callX
import static org.codehaus.groovy.ast.tools.GeneralUtils.closureX
import static org.codehaus.groovy.ast.tools.GeneralUtils.constX
import static org.codehaus.groovy.ast.tools.GeneralUtils.ctorX
import static org.codehaus.groovy.ast.tools.GeneralUtils.declX
import static org.codehaus.groovy.ast.tools.GeneralUtils.stmt
import static org.codehaus.groovy.ast.tools.GeneralUtils.varX

@CompileStatic
class ConfigAstBuilder {

    static ConfigUnit fromFileName(String filename) {
        from(CharStreams.fromFileName(filename))
    }

    static ConfigUnit fromString(String text) {
        from(CharStreams.fromString(text))
    }

    private static ConfigUnit from(CharStream inputStream) {
        final lexer = new ConfigLexer(inputStream)
        // remove the default console error listener
        // lexer.removeErrorListeners()
        final tokens = new CommonTokenStream(lexer)
        final parser = new ConfigParser(tokens)
        // remove the default console error listener
        // parser.removeErrorListeners()
        // throw an exception when a parsing error is encountered
        // parser.setErrorHandler(new BailErrorStrategy())
        parser.setBuildParseTree(true)
        final parseTree = parser.config()
        final ast = configUnit(parseTree)
        return ast
    }

    /// CONFIG

    static ConfigUnit configUnit(ConfigContext ctx) {
        final config = new ConfigUnit()
        for( final stmt : ctx.configStatement() )
            configStatement(config, stmt)
        return config
    }

    private static void configStatement(ConfigUnit config, ConfigStatementContext ctx) {
        if( ctx instanceof ConfigIncludeStmtAltContext )
            config.addStatement(configInclude(ctx.configInclude()))

        else if( ctx instanceof ConfigAssignmentStmtAltContext )
            config.addStatement(configAssignment(ctx.configAssignment()))

        else if( ctx instanceof ConfigBlockStmtAltContext )
            config.addStatement(configBlock(ctx.configBlock()))

        else
            throw new IllegalStateException()
    }

    private static Statement configInclude(ConfigIncludeContext ctx) {
        final source = expression(ctx.expression())
        final include = callX(varX('this'), 'includeConfig', argsX(source))
        stmt(include)
    }

    private static Statement configAssignment(ConfigAssignmentContext ctx) {
        final names = new ListExpression( ctx.configPathExpression().identifier().collect( ctx1 -> (Expression)constX(ctx1.text) ) )
        final right = expression(ctx.expression())
        stmt(callX(varX('this'), constX('assign'), argsX([names, right])))
    }

    private static Statement configBlock(ConfigBlockContext ctx) {
        final name = ctx.identifier()
            ? constX(ctx.identifier().text)
            : constX(unquote(ctx.stringLiteral().text))
        final statements = ctx.configBlockStatement().collect( ctx1 -> configBlockStatement(ctx1) )
        final closure = closureX(new BlockStatement(statements, new VariableScope()))
        stmt(callX(varX('this'), constX('block'), argsX([name, closure])))
    }

    private static Statement configBlockStatement(ConfigBlockStatementContext ctx) {
        if( ctx instanceof ConfigIncludeBlockStmtAltContext )
            return configInclude(ctx.configInclude())

        else if( ctx instanceof ConfigAssignmentBlockStmtAltContext )
            return configAssignment(ctx.configAssignment())

        else if( ctx instanceof ConfigBlockBlockStmtAltContext )
            return configBlock(ctx.configBlock())

        else if( ctx instanceof ConfigSelectorBlockStmtAltContext )
            return configSelector(ctx.configSelector())

        else
            throw new IllegalStateException()
    }

    private static Statement configSelector(ConfigSelectorContext ctx) {
        final kind = constX(ctx.kind.text)
        final target = configSelectorTarget(ctx.target)
        final statements = ctx.configAssignment().collect( ctx1 -> configAssignment(ctx1) )
        final closure = closureX(new BlockStatement(statements, new VariableScope()))
        stmt(callX(varX('this'), kind, argsX([target, closure])))
    }

    private static Expression configSelectorTarget(ConfigSelectorTargetContext ctx) {
        ctx.identifier()
            ? constX(ctx.identifier().text)
            : constX(unquote(ctx.stringLiteral().text))
    }

    /// STATEMENTS

    private static List<Statement> statements(List<StatementContext> stmts) {
        stmts.collect( ctx -> statement(ctx) )
    }

    static Statement statement(StatementContext ctx) {
        if( ctx instanceof BlockStmtAltContext )
            return block(ctx.block())

        if( ctx instanceof ReturnStmtAltContext )
            return returnStatement(ctx.expression())

        if( ctx instanceof AssertStmtAltContext )
            return assertStatement(ctx.assertStatement())

        if( ctx instanceof VariableDeclarationStmtAltContext ) {
            final expression = variableDeclaration(ctx.variableDeclaration())
            return stmt(expression)
        }

        if( ctx instanceof ExpressionStmtAltContext ) {
            final expression = expression(ctx.statementExpression())
            return stmt(expression)
        }

        if( ctx instanceof EmptyStmtAltContext )
            return new EmptyStatement()

        throw new IllegalStateException()
    }

    private static BlockStatement block(BlockContext ctx) {
        block(ctx.blockStatements())
    }

    private static BlockStatement block(BlockStatementsContext ctx) {
        final code = ctx
            ? statements(ctx.statement())
            : [ new EmptyStatement() as Statement ]
        new BlockStatement(code, new VariableScope())
    }

    private static ReturnStatement returnStatement(ExpressionContext ctx) {
        ctx ? new ReturnStatement(expression(ctx)) : null
    }

    private static AssertStatement assertStatement(AssertStatementContext ctx) {
        final condition = new BooleanExpression(expression(ctx.condition))
        final message = ctx.message ? expression(ctx.message) : null
        new AssertStatement(condition, message)
    }

    private static Expression variableDeclaration(VariableDeclarationContext ctx) {
        if( ctx.typeNamePairs() ) {
            final variables = ctx.typeNamePairs().typeNamePair().collect { pair ->
                final name = pair.identifier().text
                final type = type(pair.type())
                varX(name, type)
            }
            final target = variables.size() > 1
                ? new ArgumentListExpression(variables as List<Expression>)
                : variables.first()
            final init = variableInitializer(ctx.variableInitializer())
            return declX(target, init)
        }
        else {
            final type = type(ctx.type())
            final declarations = ctx.variableDeclarators().variableDeclarator().collect { decl ->
                final name = decl.identifier().text
                final target = varX(name, type)
                final init = variableInitializer(decl.variableInitializer())
                declX(target, init)
            }
            return declarations.size() > 1
                ? new ArgumentListExpression(declarations as List<Expression>)
                : declarations.first()
        }
    }

    private static Expression variableInitializer(VariableInitializerContext ctx) {
        enhancedStatementExpression(ctx.enhancedStatementExpression())
    }

    /// EXPRESSIONS

    static Expression expression(StatementExpressionContext ctx) {
        ctx.argumentList()
            ? methodCall(ctx.expression(), ctx.argumentList())
            : expression(ctx.expression())
    }

    static Expression expression(ExpressionContext ctx) {
        if( ctx instanceof AddExprAltContext )
            return binary(ctx.left, ctx.op, ctx.right)

        if( ctx instanceof AndExprAltContext )
            return binary(ctx.left, ctx.op, ctx.right)

        if( ctx instanceof AssignmentExprAltContext ) {
            final right = enhancedStatementExpression(ctx.enhancedStatementExpression())
            return binary(ctx.left, ctx.op, right)
        }

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
            final vars = ctx.variableNames().identifier().collect( ctx1 -> varX(ctx1.text) )
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
            return prefix(ctx.expression(), ctx.op)

        if( ctx instanceof UnaryNotExprAltContext )
            return prefix(ctx.expression(), ctx.op)

        throw new IllegalStateException()
    }

    private static BinaryExpression binary(ExpressionContext left, ParserToken op, ExpressionContext right) {
        new BinaryExpression(expression(left), token(op), expression(right))
    }

    private static BinaryExpression binary(ExpressionContext left, ParserToken op, Expression right) {
        new BinaryExpression(expression(left), token(op), right)
    }

    private static Expression castOperand(CastOperandExpressionContext ctx) {
        if( ctx instanceof CastCastExprAltContext ) {
            final type = type(ctx.castParExpression().type())
            final operand = castOperand(ctx.castOperandExpression())
            return new CastExpression(type, operand)
        }

        if( ctx instanceof PathCastExprAltContext )
            return path(ctx.pathExpression())

        if( ctx instanceof PostfixCastExprAltContext )
            return postfix(ctx.pathExpression(), ctx.op)

        if( ctx instanceof UnaryAddCastExprAltContext )
            return prefix(ctx.castOperandExpression(), ctx.op)

        if( ctx instanceof UnaryNotCastExprAltContext )
            return prefix(ctx.castOperandExpression(), ctx.op)

        throw new IllegalStateException()
    }

    private static Expression closure(ClosureContext ctx) {
        final params = parameters(ctx.formalParameterList())
        final code = block(ctx.blockStatements())
        final closure = closureX(params, code)
        final source = constX(ctx.text)
        ctorX(type0('ClosureWithSource'), argsX([closure, source]))
        // source
    }

    private static ConstantExpression constant(LiteralContext ctx) {
        if( ctx instanceof IntegerLiteralAltContext )
            return constX( ctx.text.toLong() )

        if( ctx instanceof FloatLiteralAltContext )
            return constX( ctx.text.toDouble() )

        if( ctx instanceof StringLiteralAltContext )
            return constX( unquote(ctx.text) )

        if( ctx instanceof BooleanLiteralAltContext )
            return constX( ctx.text=='true' ? true : false )

        if( ctx instanceof NullLiteralAltContext )
            return constX( null )

        throw new IllegalStateException()
    }

    private static Expression creator(CreatorContext ctx) {
        final type = type(ctx.createdName())
        final arguments = methodArguments(ctx.arguments().enhancedArgumentList())
        ctorX(type, arguments)
    }

    private static Expression enhancedStatementExpression(EnhancedStatementExpressionContext ctx) {
        ctx.statementExpression()
            ? expression(ctx.statementExpression())
            : lambda(ctx.standardLambdaExpression())
    }

    private static GStringExpression gstring(GstringContext ctx) {
        final strings = gstringParts(ctx)
        final values = ctx.gstringValue().collect( ctx1 -> expression(ctx1.expression()) )
        new GStringExpression(unquote(ctx.text), strings, values)
    }

    private static List<ConstantExpression> gstringParts(GstringContext ctx) {
        final strings = []

        final begin = ctx.GStringBegin().text
        strings.add(begin.startsWith('"""') ? begin[3..<-1] : begin[1..<-1])

        for( final part : ctx.GStringPart() )
            strings.add(part.text[0..<-1])

        final end = ctx.GStringEnd().text
        strings.add(begin.startsWith('"""') ? end[0..<-3] : end[0..<-1])

        return strings.collect( str -> constX(str) )
    }

    private static LambdaExpression lambda(LambdaExpressionContext ctx) {
        final params = parameters(ctx.lambdaParameters().formalParameters())
        final code = lambdaBody(ctx.lambdaBody())
        new LambdaExpression(params, code)
    }

    private static LambdaExpression lambda(StandardLambdaExpressionContext ctx) {
        final ctx1 = ctx.standardLambdaParameters()
        final params = ctx1.formalParameters()
            ? parameters(ctx1.formalParameters())
            : [ new Parameter(type(null), ctx1.identifier().text) ] as Parameter[]
        final code = lambdaBody(ctx.lambdaBody())
        new LambdaExpression(params, code)
    }

    private static BlockStatement lambdaBody(LambdaBodyContext ctx) {
        if( ctx.block() )
            return block(ctx.block())

        final statement = new ReturnStatement(expression(ctx.statementExpression()))
        new BlockStatement([ statement as Statement ], new VariableScope())
    }

    private static ListExpression list(ListContext ctx) {
        if( !ctx.expressionList() )
            return new ListExpression()
        
        final expressions = ctx.expressionList().expressionListElement().collect( this.&listElement )
        new ListExpression(expressions)
    }

    private static Expression listElement(ExpressionListElementContext ctx) {
        final element = expression(ctx.expression())
        ctx.MUL()
            ? new SpreadExpression(element)
            : element
    }

    private static MapExpression map(MapContext ctx) {
        if( !ctx.mapEntryList() )
            return new MapExpression()

        final entries = ctx.mapEntryList().mapEntry().collect( this.&mapEntry )
        new MapExpression(entries)
    }

    private static MapEntryExpression mapEntry(MapEntryContext ctx) {
        final value = expression(ctx.expression())
        final key = ctx.MUL()
            ? new SpreadMapExpression(value)
            : mapEntryLabel(ctx.mapEntryLabel())
        new MapEntryExpression(key, value)
    }

    private static Expression mapEntryLabel(MapEntryLabelContext ctx) {
        if( ctx.keywords() )
            return constX(ctx.text)

        if( ctx.primary() )
            return primary(ctx.primary())

        throw new IllegalStateException()
    }

    private static MethodCallExpression methodCall(ExpressionContext method, ArgumentListContext args) {
        final parts = methodObject(expression(method))
        final arguments = methodArguments(args)
        callX(parts[0], parts[1], arguments)
    }

    private static List<Expression> methodObject(Expression expression) {
        if( expression instanceof PropertyExpression )
            return [expression.objectExpression, expression.property]

        if( expression instanceof VariableExpression )
            return [varX('this'), constX(expression.text)]

        return [expression, constX('call')]
    }

    private static Expression methodArguments(ArgumentListContext ctx) {
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
                throw new IllegalStateException()
        }

        if( opts )
            args.push(new MapExpression(opts))

        return new ArgumentListExpression(args)
    }

    private static Expression methodArguments(EnhancedArgumentListContext ctx) {
        if( !ctx )
            return new ArgumentListExpression()

        final List<Expression> args = []
        final List<MapEntryExpression> opts = []

        for( final ctx1 : ctx.enhancedArgumentListElement() ) {
            if( ctx1.expressionListElement() )
                args << listElement(ctx1.expressionListElement())

            else if( ctx1.standardLambdaExpression() )
                args << lambda(ctx1.standardLambdaExpression())

            else if( ctx1.namedArg() )
                opts << namedArg(ctx1.namedArg())

            else
                throw new IllegalStateException()
        }

        if( opts )
            args.push(new MapExpression(opts))

        return new ArgumentListExpression(args)
    }

    private static MapEntryExpression namedArg(NamedArgContext ctx) {
        final value = expression(ctx.expression())
        final key = ctx.MUL()
            ? new SpreadMapExpression(value)
            : namedArgLabel(ctx.namedArgLabel())
        new MapEntryExpression(key, value)
    }

    private static Expression namedArgLabel(NamedArgLabelContext ctx) {
        if( ctx.keywords() || ctx.identifier() )
            return constX(ctx.text)

        if( ctx.literal() )
            return constant(ctx.literal())

        if( ctx.gstring() )
            return gstring(ctx.gstring())

        throw new IllegalStateException()
    }

    private static Expression path(PathExpressionContext ctx) {
        def result = primary(ctx.primary())
        for( final el : ctx.pathElement() )
            result = pathElement(result, el)
        return result
    }

    private static Expression pathElement(Expression expression, PathElementContext ctx) {
        if( ctx instanceof PropertyPathExprAltContext ) {
            final prop =
                ctx.keywords() ? ctx.keywords().text :
                ctx.identifier() ? ctx.identifier().text :
                /* ctx.stringLiteral() */ unquote(ctx.stringLiteral().text)
            final safe = ctx.SAFE_DOT() != null || ctx.SPREAD_DOT() != null
            final result = ctx.AT()
                ? new AttributeExpression(expression, constX(prop), safe)
                : new PropertyExpression(expression, constX(prop), safe)
            if( ctx.SPREAD_DOT() )
                result.setSpreadSafe(true)
            return result
        }

        if( ctx instanceof MethodCallPathExprAltContext ) {
            final parts = methodObject(expression)
            final arguments = methodArguments(ctx.arguments().enhancedArgumentList())
            return callX(parts[0], parts[1], arguments)
        }

        if( ctx instanceof ListElementPathExprAltContext ) {
            final op = new Token(Types.LEFT_SQUARE_BRACKET, '[', -1, -1)
            final elements = ctx.expressionList().collect( this.&expression )
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

    private static PostfixExpression postfix(PathExpressionContext ctx, ParserToken op) {
        new PostfixExpression(path(ctx), token(op))
    }

    private static PrefixExpression prefix(ExpressionContext ctx, ParserToken op) {
        new PrefixExpression(token(op), expression(ctx))
    }

    private static PrefixExpression prefix(CastOperandExpressionContext ctx, ParserToken op) {
        new PrefixExpression(token(op), castOperand(ctx))
    }

    private static Expression primary(PrimaryContext ctx) {
        if( ctx instanceof IdentifierPrmrAltContext )
            return varX(ctx.identifier().text)

        if( ctx instanceof LiteralPrmrAltContext )
            return constant(ctx.literal())

        if( ctx instanceof GstringPrmrAltContext )
            return gstring(ctx.gstring())

        if( ctx instanceof NewPrmrAltContext )
            return creator(ctx.creator())

        if( ctx instanceof ParenPrmrAltContext )
            return enhancedStatementExpression(ctx.parExpression().enhancedStatementExpression())

        if( ctx instanceof ClosureOrLambdaPrmrAltContext ) {
            final ctx1 = ctx.closureOrLambdaExpression()
            return ctx1.closure()
                ? closure(ctx1.closure())
                : lambda(ctx1.lambdaExpression())
        }

        if( ctx instanceof ListPrmrAltContext )
            return list(ctx.list())

        if( ctx instanceof MapPrmrAltContext )
            return map(ctx.map())

        if( ctx instanceof BuiltInTypePrmrAltContext )
            return new ClassExpression(type(ctx.builtInType()))

        throw new IllegalStateException()
    }

    private static TernaryExpression ternary(ExpressionContext cond, ExpressionContext tb, ExpressionContext fb) {
        final condition = new BooleanExpression(expression(cond))
        final trueExpression = tb ? expression(tb) : null
        final falseExpression = expression(fb)
        new TernaryExpression(condition, trueExpression, falseExpression)
    }

    /// MISCELLANEOUS

    private static Parameter[] parameters(FormalParametersContext ctx) {
        parameters(ctx.formalParameterList())
    }

    private static Parameter[] parameters(FormalParameterListContext ctx) {
        final result = ctx
            ? ctx.formalParameter().collect( this.&parameter )
            : []
        return result as Parameter[]
    }

    private static Parameter parameter(FormalParameterContext ctx) {
        final type = type(ctx.type())
        final name = ctx.identifier().text
        final defaultValue = ctx.expression() ? expression(ctx.expression()) : null
        new Parameter(type, name, defaultValue)
    }

    private static Token token(ParserToken t) {
        new Token(Types.lookupSymbol(t.text), t.text, -1, -1)
    }

    private static ClassNode type(ParserRuleContext ctx) {
        ctx
            ? type0(ctx.text)
            : ClassHelper.DYNAMIC_TYPE
    }

    private static ClassNode type0(String text) {
        ClassHelper.makeWithoutCaching(text)
    }

    /// HELPERS

    private static String unquote(String text) {
        if( !text )
            return null
        else if( text.startsWith('"""') || text.startsWith("'''") )
            return text[3..<-3]
        else
            return text[1..<-1]
    }

}

@CompileStatic
class ConfigUnit {
    private List<Statement> statements = []

    void addStatement(Statement stmt) {
        statements.add(stmt)
    }

    List<Statement> getStatements() {
        statements
    }
}
