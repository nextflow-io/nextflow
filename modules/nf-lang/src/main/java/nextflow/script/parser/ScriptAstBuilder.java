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
package nextflow.script.parser;

import java.io.BufferedReader;
import java.io.IOException;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.atomic.AtomicInteger;

import groovy.lang.Tuple2;
import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.ast.AssignmentExpression;
import nextflow.script.ast.FeatureFlagNode;
import nextflow.script.ast.FunctionNode;
import nextflow.script.ast.IncludeModuleNode;
import nextflow.script.ast.IncludeNode;
import nextflow.script.ast.IncompleteNode;
import nextflow.script.ast.InvalidDeclaration;
import nextflow.script.ast.OutputNode;
import nextflow.script.ast.ParamNode;
import nextflow.script.ast.ProcessNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.WorkflowNode;
import org.antlr.v4.runtime.ANTLRErrorListener;
import org.antlr.v4.runtime.CharStream;
import org.antlr.v4.runtime.CharStreams;
import org.antlr.v4.runtime.CommonTokenStream;
import org.antlr.v4.runtime.Parser;
import org.antlr.v4.runtime.ParserRuleContext;
import org.antlr.v4.runtime.RecognitionException;
import org.antlr.v4.runtime.Recognizer;
import org.antlr.v4.runtime.Token;
import org.antlr.v4.runtime.atn.ATNConfigSet;
import org.antlr.v4.runtime.atn.PredictionMode;
import org.antlr.v4.runtime.dfa.DFA;
import org.antlr.v4.runtime.misc.ParseCancellationException;
import org.antlr.v4.runtime.tree.TerminalNode;
import org.apache.groovy.parser.antlr4.GroovySyntaxError;
import org.apache.groovy.parser.antlr4.util.StringUtils;
import org.codehaus.groovy.GroovyBugError;
import org.codehaus.groovy.antlr.EnumHelper;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.GenericsType;
import org.codehaus.groovy.ast.ModuleNode;
import org.codehaus.groovy.ast.NodeMetaDataHandler;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.VariableScope;
import org.codehaus.groovy.ast.expr.ArgumentListExpression;
import org.codehaus.groovy.ast.expr.BinaryExpression;
import org.codehaus.groovy.ast.expr.BitwiseNegationExpression;
import org.codehaus.groovy.ast.expr.ClassExpression;
import org.codehaus.groovy.ast.expr.ConstantExpression;
import org.codehaus.groovy.ast.expr.EmptyExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.GStringExpression;
import org.codehaus.groovy.ast.expr.MapExpression;
import org.codehaus.groovy.ast.expr.MapEntryExpression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.NamedArgumentListExpression;
import org.codehaus.groovy.ast.expr.NotExpression;
import org.codehaus.groovy.ast.expr.PropertyExpression;
import org.codehaus.groovy.ast.expr.RangeExpression;
import org.codehaus.groovy.ast.expr.TupleExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.expr.UnaryMinusExpression;
import org.codehaus.groovy.ast.expr.UnaryPlusExpression;
import org.codehaus.groovy.ast.stmt.AssertStatement;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.CatchStatement;
import org.codehaus.groovy.ast.stmt.EmptyStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.CompilationFailedException;
import org.codehaus.groovy.control.CompilePhase;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;
import org.codehaus.groovy.syntax.Numbers;
import org.codehaus.groovy.syntax.SyntaxException;
import org.codehaus.groovy.syntax.Types;

import static nextflow.script.ast.ASTUtils.*;
import static nextflow.script.parser.PositionConfigureUtils.ast;
import static nextflow.script.parser.PositionConfigureUtils.tokenPosition;
import static nextflow.script.parser.ScriptParser.*;
import static org.codehaus.groovy.ast.expr.VariableExpression.THIS_EXPRESSION;
import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Transform a Nextflow script parse tree into an abstract syntax tree (AST).
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ScriptAstBuilder {

    private SourceUnit sourceUnit;
    private ScriptNode moduleNode;
    private ScriptLexer lexer;
    private ScriptParser parser;
    private final GroovydocManager groovydocManager;

    private Tuple2<ParserRuleContext,Exception> numberFormatError;

    public ScriptAstBuilder(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
        this.moduleNode = new ScriptNode(sourceUnit);

        var charStream = createCharStream(sourceUnit);
        this.lexer = new ScriptLexer(charStream);
        this.parser = new ScriptParser(new CommonTokenStream(lexer));
        parser.setErrorHandler(new DescriptiveErrorStrategy(charStream));

        var groovydocEnabled = sourceUnit.getConfiguration().isGroovydocEnabled();
        this.groovydocManager = new GroovydocManager(groovydocEnabled);
    }

    private CharStream createCharStream(SourceUnit sourceUnit) {
        try {
            return CharStreams.fromReader(
                    new BufferedReader(sourceUnit.getSource().getReader()),
                    sourceUnit.getName());
        }
        catch( IOException e ) {
            throw new RuntimeException("Error occurred when reading source code.", e);
        }
    }

    private CompilationUnitContext buildCST() {
        try {
            var tokenStream = parser.getInputStream();
            try {
                return buildCST(PredictionMode.SLL);
            }
            catch( Throwable t ) {
                // if some syntax error occurred in the lexer, no need to retry the powerful LL mode
                if( t instanceof GroovySyntaxError gse && gse.getSource() == GroovySyntaxError.LEXER )
                    throw t;

                // System.err.println("SLL parsing failed: " + sourceUnit.getSource().getURI().toString());
                tokenStream.seek(0);
                return buildCST(PredictionMode.LL);
            }
        }
        catch( Throwable t ) {
            throw convertException(t);
        }
    }

    private CompilationUnitContext buildCST(PredictionMode predictionMode) {
        parser.getInterpreter().setPredictionMode(predictionMode);

        removeErrorListeners();
        if( predictionMode == PredictionMode.LL )
            addErrorListeners();

        return parser.compilationUnit();
    }

    private CompilationFailedException convertException(Throwable t) {
        if( t instanceof CompilationFailedException cfe )
            return cfe;
        else if( t instanceof ParseCancellationException )
            return createParsingFailedException(t.getCause());
        else
            return createParsingFailedException(t);
    }

    public ModuleNode buildAST() {
        try {
            return compilationUnit(buildCST());
        }
        catch( Throwable t ) {
            throw convertException(t);
        }
    }

    /// SCRIPT DECLARATIONS

    private static final List<String> SCRIPT_DEF_NAMES = List.of("process", "workflow", "output");

    private ModuleNode compilationUnit(CompilationUnitContext ctx) {
        var statements = new ArrayList<Statement>();
        boolean hasDeclarations = false;

        for( var declOrStmt : ctx.scriptDeclarationOrStatement() ) {
            if( declOrStmt.scriptDeclaration() != null )
                hasDeclarations |= scriptDeclaration(declOrStmt.scriptDeclaration());
            if( declOrStmt.statement() != null )
                statements.add(statement(declOrStmt.statement()));
        }

        for( int i = 0; i < statements.size(); i++ ) {
            var stmt = statements.get(i);
            var defName = stmt instanceof ExpressionStatement es && es.getExpression() instanceof MethodCallExpression mce
                ? mce.getMethodAsString()
                : null;
            if( defName != null && SCRIPT_DEF_NAMES.contains(defName) ) {
                collectSyntaxError(new SyntaxException("Invalid " + defName + " definition -- check for missing or out-of-order section labels", stmt));
                statements.set(i, ast(new InvalidDeclaration(), stmt));
                hasDeclarations = true;
            }
        }

        if( hasDeclarations ) {
            for( var stmt : statements ) {
                if( !(stmt instanceof InvalidDeclaration) )
                    collectSyntaxError(new SyntaxException("Statements cannot be mixed with script declarations -- move statements into a process or workflow", stmt));
            }
        }

        if( !statements.isEmpty() ) {
            var main = block(new VariableScope(), statements);
            ast( main, statements.get(0), statements.get(statements.size() - 1) );
            var workflowNode = workflowDef(main);
            moduleNode.addWorkflow(workflowNode);
            if( !hasDeclarations )
                moduleNode.setEntry(workflowNode);
        }

        // NOTE: required to produce a valid script class
        moduleNode.addStatement(EmptyStatement.INSTANCE);

        if( numberFormatError != null )
            throw createParsingFailedException(numberFormatError.getV2().getMessage(), numberFormatError.getV1());

        return moduleNode;
    }

    private boolean scriptDeclaration(ScriptDeclarationContext ctx) {
        if( ctx instanceof FeatureFlagDeclAltContext ffac ) {
            var node = featureFlagDeclaration(ffac.featureFlagDeclaration());
            saveLeadingComments(node, ctx);
            moduleNode.addFeatureFlag(node);
        }

        else if( ctx instanceof EnumDefAltContext edac ) {
            var node = enumDef(edac.enumDef());
            saveLeadingComments(node, ctx);
            moduleNode.addClass(node);
        }

        else if( ctx instanceof FunctionDefAltContext fdac ) {
            var node = functionDef(fdac.functionDef());
            saveLeadingComments(node, ctx);
            moduleNode.addFunction(node);
        }

        else if( ctx instanceof ImportDeclAltContext iac ) {
            var node = ast( new EmptyStatement(), iac );
            collectSyntaxError(new SyntaxException("Groovy `import` declarations are not supported -- use fully-qualified name inline instead", node));
            return false;
        }

        else if( ctx instanceof IncludeDeclAltContext iac ) {
            var node = includeDeclaration(iac.includeDeclaration());
            saveLeadingComments(node, ctx);
            moduleNode.addInclude(node);
        }

        else if( ctx instanceof OutputDefAltContext odac ) {
            var node = outputDef(odac.outputDef());
            saveLeadingComments(node, ctx);
            if( moduleNode.getOutput() != null )
                collectSyntaxError(new SyntaxException("Output block defined more than once", node));
            moduleNode.setOutput(node);
        }

        else if( ctx instanceof ParamDeclAltContext pac ) {
            var node = paramDeclaration(pac.paramDeclaration());
            saveLeadingComments(node, ctx);
            moduleNode.addParam(node);
        }

        else if( ctx instanceof ProcessDefAltContext pdac ) {
            var node = processDef(pdac.processDef());
            saveLeadingComments(node, ctx);
            moduleNode.addProcess(node);
        }

        else if( ctx instanceof WorkflowDefAltContext wdac ) {
            var node = workflowDef(wdac.workflowDef());
            saveLeadingComments(node, ctx);
            if( node.getName() == null ) {
                if( moduleNode.getEntry() != null )
                    collectSyntaxError(new SyntaxException("Entry workflow defined more than once", node));
                moduleNode.setEntry(node);
            }
            moduleNode.addWorkflow(node);
        }

        else if( ctx instanceof IncompleteScriptDeclAltContext iac ) {
            incompleteScriptDeclaration(iac.incompleteScriptDeclaration());
            return false;
        }

        else
            throw createParsingFailedException("Invalid script declaration: " + ctx.getText(), ctx);

        return true;
    }

    private FeatureFlagNode featureFlagDeclaration(FeatureFlagDeclarationContext ctx) {
        var name = ctx.featureFlagName().getText();
        var value = expression(ctx.expression());
        var result = ast( new FeatureFlagNode(name, value), ctx );
        if( !(value instanceof ConstantExpression) )
            collectSyntaxError(new SyntaxException("Feature flag value must be a literal value (number, string, true/false)", result));
        return result;
    }

    private ParamNode paramDeclaration(ParamDeclarationContext ctx) {
        Expression target = ast( varX("params"), ctx.PARAMS() );
        for( var ident : ctx.identifier() ) {
            var name = ast( constX(identifier(ident)), ident );
            target = ast( propX(target, name), target, name );
        }
        var value = expression(ctx.expression());
        return ast( new ParamNode(target, value), ctx );
    }

    private IncludeNode includeDeclaration(IncludeDeclarationContext ctx) {
        var source = ast( string(ctx.stringLiteral()), ctx.stringLiteral() );
        var modules = ctx.includeNames().includeName().stream()
            .map(it -> {
                var name = it.name.getText();
                var alias = it.alias != null ? it.alias.getText() : null;
                var result = new IncludeModuleNode(name, alias);
                result.putNodeMetaData("_START_NAME", tokenPosition(it.name));
                if( it.alias != null )
                    result.putNodeMetaData("_START_ALIAS", tokenPosition(it.alias));
                return ast( result, it );
            })
            .toList();

        return ast( new IncludeNode(source, modules), ctx );
    }

    private ClassNode enumDef(EnumDefContext ctx) {
        var name = identifier(ctx.identifier());
        var result = ast( EnumHelper.makeEnumNode(name, Modifier.PUBLIC, null, null), ctx );
        if( ctx.enumBody() != null ) {
            for( var ident : ctx.enumBody().identifier() ) {
                var fn = ast( EnumHelper.addEnumConstant(result, identifier(ident), null), ident );
                groovydocManager.handle(fn, ident);
            }
        }
        groovydocManager.handle(result, ctx);
        return result;
    }

    private ProcessNode processDef(ProcessDefContext ctx) {
        var name = ctx.name.getText();
        if( ctx.body == null ) {
            var empty = EmptyStatement.INSTANCE;
            var result = ast( new ProcessNode(name, empty, empty, empty, EmptyExpression.INSTANCE, null, empty, empty), ctx );
            collectSyntaxError(new SyntaxException("Missing process script body", result));
            return result;
        }

        var directives = processDirectives(ctx.body.processDirectives());
        var inputs = processInputs(ctx.body.processInputs());
        var outputs = processOutputs(ctx.body.processOutputs());
        var when = processWhen(ctx.body.processWhen());
        var type = processType(ctx.body.processExec());
        var exec = ctx.body.blockStatements() != null
            ? blockStatements(ctx.body.blockStatements())
            : blockStatements(ctx.body.processExec().blockStatements());
        var stub = processStub(ctx.body.processStub());

        if( ctx.body.blockStatements() != null ) {
            if( !(directives instanceof EmptyStatement) || !(inputs instanceof EmptyStatement) || !(outputs instanceof EmptyStatement) )
                collectSyntaxError(new SyntaxException("The `script:` or `exec:` label is required when other sections are present", exec));
        }

        var result = ast( new ProcessNode(name, directives, inputs, outputs, when, type, exec, stub), ctx );
        groovydocManager.handle(result, ctx);
        return result;
    }

    private Statement processDirectives(ProcessDirectivesContext ctx) {
        if( ctx == null )
            return EmptyStatement.INSTANCE;
        var statements = ctx.statement().stream()
            .map(this::statement)
            .map(stmt -> checkDirective(stmt, "Invalid process directive"))
            .toList();
        return ast( block(null, statements), ctx );
    }

    private Statement processInputs(ProcessInputsContext ctx) {
        if( ctx == null )
            return EmptyStatement.INSTANCE;
        var statements = ctx.statement().stream()
            .map(this::statement)
            .map(stmt -> checkDirective(stmt, "Invalid process input"))
            .toList();
        return ast( block(null, statements), ctx );
    }

    private Statement processOutputs(ProcessOutputsContext ctx) {
        if( ctx == null )
            return EmptyStatement.INSTANCE;
        var statements = ctx.statement().stream()
            .map(this::statement)
            .map(stmt -> checkDirective(stmt, "Invalid process output"))
            .toList();
        return ast( block(null, statements), ctx );
    }

    private Statement checkDirective(Statement stmt, String errorMessage) {
        if( !(stmt instanceof ExpressionStatement) ) {
            collectSyntaxError(new SyntaxException(errorMessage, stmt));
            return ast( new EmptyStatement(), stmt );
        }

        var stmtX = (ExpressionStatement) stmt;
        var expression = stmtX.getExpression();
        if( expression instanceof VariableExpression ve ) {
            var method = ast( constX(ve.getName()), ve );
            var call = ast( callX(THIS_EXPRESSION, method, new ArgumentListExpression()), stmtX );
            stmtX.setExpression(call);
            return stmtX;
        }
        if( isDirectiveWithNegativeValue(expression) ) {
            var binary = (BinaryExpression) expression;
            var left = (VariableExpression) binary.getLeftExpression();
            var method = ast( constX(left.getName()), left );
            var value = (Expression) ast( new UnaryMinusExpression(binary.getRightExpression()), binary.getRightExpression() );
            var arguments = ast( args(value), value );
            var call = ast( callX(THIS_EXPRESSION, method, arguments), stmtX );
            stmtX.setExpression(call);
            return stmtX;
        }
        if( !(expression instanceof MethodCallExpression) ) {
            collectSyntaxError(new SyntaxException(errorMessage, stmtX));
            return ast( new EmptyStatement(), stmtX );
        }

        var call = (MethodCallExpression) expression;
        if( !call.isImplicitThis() || !(call.getMethod() instanceof ConstantExpression) ) {
            collectSyntaxError(new SyntaxException(errorMessage, stmtX));
            return ast( new EmptyStatement(), stmtX );
        }

        return stmtX;
    }

    private boolean isDirectiveWithNegativeValue(Expression expression) {
        if( !(expression instanceof BinaryExpression) )
            return false;
        var binary = (BinaryExpression) expression;
        if( !(binary.getLeftExpression() instanceof VariableExpression) )
            return false;
        if( binary.getOperation().getType() != Types.MINUS )
            return false;
        return true;
    }

    private Expression processWhen(ProcessWhenContext ctx) {
        if( ctx == null )
            return EmptyExpression.INSTANCE;
        return ast( expression(ctx.expression()), ctx );
    }

    private String processType(ProcessExecContext ctx) {
        if( ctx == null )
            return "script";

        if( ctx.EXEC() != null ) {
            return "exec";
        }
        if( ctx.SHELL() != null ) {
            collectWarning("The `shell` block is deprecated, use `script` instead", ast( new EmptyExpression(), ctx.SHELL() ));
            return "shell";
        }
        return "script";
    }

    private Statement processStub(ProcessStubContext ctx) {
        if( ctx == null )
            return EmptyStatement.INSTANCE;
        return ast( blockStatements(ctx.blockStatements()), ctx );
    }

    private WorkflowNode workflowDef(WorkflowDefContext ctx) {
        var name = ctx.name != null ? ctx.name.getText() : null;

        if( ctx.body == null ) {
            var result = ast( new WorkflowNode(name, null, null, null, null), ctx );
            groovydocManager.handle(result, ctx);
            return result;
        }

        var takes = workflowTakes(ctx.body.workflowTakes());
        var emits = workflowEmits(ctx.body.workflowEmits());
        var publishers = workflowPublishers(ctx.body.workflowPublishers());
        var main = blockStatements(
            ctx.body.workflowMain() != null
                ? ctx.body.workflowMain().blockStatements()
                : null
        );

        if( name == null ) {
            if( takes instanceof BlockStatement )
                collectSyntaxError(new SyntaxException("Entry workflow cannot have a take section", takes));
            if( emits instanceof BlockStatement )
                collectSyntaxError(new SyntaxException("Entry workflow cannot have an emit section", emits));
        }

        var result = ast( new WorkflowNode(name, takes, main, emits, publishers), ctx );
        groovydocManager.handle(result, ctx);
        return result;
    }

    private WorkflowNode workflowDef(BlockStatement main) {
        var takes = EmptyStatement.INSTANCE;
        var emits = EmptyStatement.INSTANCE;
        var publishers = EmptyStatement.INSTANCE;
        return new WorkflowNode(null, takes, main, emits, publishers);
    }

    private Statement workflowTakes(WorkflowTakesContext ctx) {
        if( ctx == null )
            return EmptyStatement.INSTANCE;

        var statements = ctx.identifier().stream()
            .map(this::workflowTake)
            .toList();
        return ast( block(null, statements), ctx );
    }

    private Statement workflowTake(IdentifierContext ctx) {
        var result = ast( stmt(variableName(ctx)), ctx );
        saveTrailingComment(result, ctx);
        return result;
    }

    private Statement workflowEmits(WorkflowEmitsContext ctx) {
        if( ctx == null )
            return EmptyStatement.INSTANCE;

        var statements = ctx.statement().stream()
            .map(this::workflowEmit)
            .filter(stmt -> stmt != null)
            .toList();
        var result = ast( block(null, statements), ctx );
        var hasEmitExpression = statements.stream().anyMatch(this::isEmitExpression);
        if( hasEmitExpression && statements.size() > 1 )
            collectSyntaxError(new SyntaxException("Every emit must be assigned to a name when there are multiple emits", result));
        return result;
    }

    private Statement workflowEmit(StatementContext ctx) {
        var result = statement(ctx);
        if( !(result instanceof ExpressionStatement) ) {
            collectSyntaxError(new SyntaxException("Invalid workflow emit -- must be a name, assignment, or expression", result));
            return null;
        }
        saveTrailingComment(result, ctx);
        return result;
    }

    private boolean isEmitExpression(Statement stmt) {
        if( stmt instanceof ExpressionStatement es ) {
            var exp = es.getExpression();
            return !(exp instanceof VariableExpression || exp instanceof AssignmentExpression);
        }
        return false;
    }

    private Statement workflowPublishers(WorkflowPublishersContext ctx) {
        if( ctx == null )
            return EmptyStatement.INSTANCE;

        var statements = ctx.statement().stream()
            .map(this::statement)
            .map(this::checkWorkflowPublisher)
            .filter(stmt -> stmt != null)
            .toList();
        return ast( block(null, statements), ctx );
    }

    private Statement checkWorkflowPublisher(Statement stmt) {
        var valid = stmt instanceof ExpressionStatement es
            && es.getExpression() instanceof BinaryExpression be
            && be.getOperation().getType() == Types.RIGHT_SHIFT;
        if( !valid ) {
            collectSyntaxError(new SyntaxException("Invalid workflow publish statement", stmt));
            return null;
        }
        return stmt;
    }

    private OutputNode outputDef(OutputDefContext ctx) {
        var body = outputBody(ctx.outputBody());
        return ast( new OutputNode(body), ctx );
    }

    private Statement outputBody(OutputBodyContext ctx) {
        if( ctx == null )
            return EmptyStatement.INSTANCE;
        var statements = ctx.outputTargetBody().stream()
            .map(this::outputTargetBody)
            .toList();
        return ast( block(null, statements), ctx );
    }

    private Statement outputTargetBody(OutputTargetBodyContext ctx) {
        var stmt = statement(ctx.statement());
        var call = asMethodCallX(stmt);
        if( call != null ) {
            var block = asDslBlock(call, 1);
            if( block != null )
                return stmt;
        }
        collectSyntaxError(new SyntaxException("Invalid output target block", stmt));
        return ast( new EmptyStatement(), stmt );
    }

    private FunctionNode functionDef(FunctionDefContext ctx) {
        var name = identifier(ctx.identifier());
        var returnType = legacyType(ctx.legacyType());
        var params = Optional.ofNullable(formalParameterList(ctx.formalParameterList())).orElse(Parameter.EMPTY_ARRAY);
        var code = blockStatements(ctx.blockStatements());

        var result = ast( new FunctionNode(name, returnType, params, code), ctx );
        checkInvalidVarName(name, result);
        groovydocManager.handle(result, ctx);
        return result;
    }

    private Statement incompleteScriptDeclaration(IncompleteScriptDeclarationContext ctx) {
        var result = ast( new IncompleteNode(ctx.getText()), ctx );
        collectSyntaxError(new SyntaxException("Incomplete declaration", result));
        return result;
    }

    /// STATEMENTS

    private Statement statement(StatementContext ctx) {
        Statement result;

        if( ctx instanceof IfElseStmtAltContext ieac )
            result = ast( ifElseStatement(ieac.ifElseStatement()), ieac );

        else if( ctx instanceof TryCatchStmtAltContext tcac )
            result = ast( tryCatchStatement(tcac.tryCatchStatement()), tcac );

        else if( ctx instanceof ReturnStmtAltContext rac )
            result = ast( returnStatement(rac.expression()), rac );

        else if( ctx instanceof ThrowStmtAltContext tac )
            result = ast( throwStatement(tac.expression()), tac );

        else if( ctx instanceof AssertStmtAltContext aac )
            result = ast( assertStatement(aac.assertStatement()), aac );

        else if( ctx instanceof VariableDeclarationStmtAltContext vdac )
            result = ast( variableDeclaration(vdac.variableDeclaration()), vdac );

        else if( ctx instanceof MultipleAssignmentStmtAltContext maac )
            result = ast( assignment(maac.multipleAssignmentStatement()), maac );

        else if( ctx instanceof AssignmentStmtAltContext aac )
            result = ast( assignment(aac.assignmentStatement()), aac );

        else if( ctx instanceof ExpressionStmtAltContext eac )
            result = ast( expressionStatement(eac.expressionStatement()), eac );

        else if( ctx instanceof EmptyStmtAltContext )
            return EmptyStatement.INSTANCE;

        else
            throw createParsingFailedException("Invalid statement: " + ctx.getText(), ctx);

        saveLeadingComments(result, ctx);
        return result;
    }

    private Statement ifElseStatement(IfElseStatementContext ctx) {
        var expression = ast( parExpression(ctx.parExpression()), ctx.parExpression() );
        var condition = ast( boolX(expression), expression );
        var thenStmt = statementOrBlock(ctx.tb);
        var elseStmt = ctx.ELSE() != null
            ? statementOrBlock(ctx.fb)
            : EmptyStatement.INSTANCE;
        return ifElseS(condition, thenStmt, elseStmt);
    }

    private Statement statementOrBlock(StatementOrBlockContext ctx) {
        return ctx.statement() != null
            ? statement(ctx.statement())
            : blockStatements(ctx.blockStatements());
    }

    private BlockStatement blockStatements(BlockStatementsContext ctx) {
        if( ctx == null )
            return block(new VariableScope(), Collections.emptyList());
        var statements = ctx.statement().stream()
            .map(this::statement)
            .toList();
        return ast( block(new VariableScope(), statements), ctx );
    }

    private Statement tryCatchStatement(TryCatchStatementContext ctx) {
        var tryStatement = statementOrBlock(ctx.statementOrBlock());
        var catchClauses = ctx.catchClause().stream()
            .map(this::catchClause)
            .toList();
        var result = tryCatchS(tryStatement);
        for( var clause : catchClauses )
            for( var stmt : clause )
                result.addCatch(stmt);
        return result;
    }

    private List<CatchStatement> catchClause(CatchClauseContext ctx) {
        var types = catchTypes(ctx.catchTypes());
        return types.stream()
            .map(type -> {
                var name = identifier(ctx.identifier());
                var variable = ast( param(type, name), ctx.identifier() );
                checkInvalidVarName(name, variable);
                var code = statementOrBlock(ctx.statementOrBlock());
                return ast( new CatchStatement(variable, code), ctx );
            })
            .toList();
    }

    private List<ClassNode> catchTypes(CatchTypesContext ctx) {
        if( ctx == null )
            return Collections.singletonList( ClassHelper.dynamicType() );

        return ctx.qualifiedClassName().stream()
            .map(this::qualifiedClassName)
            .toList();
    }

    private Statement returnStatement(ExpressionContext ctx) {
        var result = ctx != null
            ? expression(ctx)
            : ConstantExpression.EMPTY_EXPRESSION;
        return returnS(result);
    }

    private Statement throwStatement(ExpressionContext ctx) {
        var result = ctx != null
            ? expression(ctx)
            : ConstantExpression.EMPTY_EXPRESSION;
        return throwS(result);
    }

    private Statement assertStatement(AssertStatementContext ctx) {
        var condition = ast( boolX(expression(ctx.condition)), ctx.condition );
        return ctx.message != null
            ? new AssertStatement(condition, expression(ctx.message))
            : new AssertStatement(condition);
    }

    private Statement variableDeclaration(VariableDeclarationContext ctx) {
        if( ctx.variableNames() != null ) {
            // multiple assignment
            var variables = ctx.variableNames().identifier().stream()
                .map(ident -> (Expression) variableName(ident))
                .toList();
            var target = new ArgumentListExpression(variables);
            var initializer = expression(ctx.initializer);
            return stmt(ast( declX(target, initializer), ctx ));
        }
        else {
            // single assignment
            var target = variableName(ctx.identifier());
            target.setType(legacyType(ctx.legacyType()));
            var initializer = ctx.initializer != null
                ? expression(ctx.initializer)
                : EmptyExpression.INSTANCE;
            return stmt(ast( declX(target, initializer), ctx ));
        }
    }

    private Expression variableNames(VariableNamesContext ctx) {
        var vars = ctx.identifier().stream()
            .map(this::variableName)
            .toList();
        return ast( new TupleExpression(vars), ctx );
    }

    private Expression variableName(IdentifierContext ctx) {
        var name = identifier(ctx);
        var result = ast( varX(name), ctx );
        checkInvalidVarName(name, result);
        return result;
    }

    private static final String KEYWORD_FOR = Types.getText(Types.KEYWORD_FOR);
    private static final String KEYWORD_SWITCH = Types.getText(Types.KEYWORD_SWITCH);
    private static final String KEYWORD_WHILE = Types.getText(Types.KEYWORD_WHILE);

    private void checkInvalidVarName(String name, ASTNode node) {
        if( "_".equals(name) )
            collectSyntaxError(new SyntaxException("`_` is not allowed as an identifier because it is reserved for future use", node));
        if( !GROOVY_KEYWORDS.contains(name) )
            return;
        if( KEYWORD_FOR.equals(name) || KEYWORD_WHILE.equals(name) )
            collectSyntaxError(new SyntaxException("`" + name + "` loops are no longer supported", node));
        else if( KEYWORD_SWITCH.equals(name) )
            collectSyntaxError(new SyntaxException("`switch` statements are no longer supported", node));
        else
            collectSyntaxError(new SyntaxException("`" + name + "` is not allowed as an identifier because it is a Groovy keyword", node));
    }

    private Statement assignment(MultipleAssignmentStatementContext ctx) {
        var target = variableNames(ctx.variableNames());
        var source = expression(ctx.expression());
        return stmt(ast( new AssignmentExpression(target, source), ctx ));
    }

    private Statement assignment(AssignmentStatementContext ctx) {
        var target = expression(ctx.target);
        if( target instanceof VariableExpression && isInsideParentheses(target) ) {
            if( target.<Number>getNodeMetaData(ASTNodeMarker.INSIDE_PARENTHESES_LEVEL).intValue() > 1 )
                throw createParsingFailedException("Nested parenthesis is not allowed in multiple assignment, e.g. ((a)) = b", ctx);

            var tuple = ast( new TupleExpression(target), ctx.target );
            return stmt(ast( new AssignmentExpression(tuple, token(ctx.op), expression(ctx.source)), ctx ));
        }

        if ( isValidAssignmentTarget(target) )
            return stmt(ast( new AssignmentExpression(target, token(ctx.op), expression(ctx.source)), ctx ));

        throw createParsingFailedException("Invalid assignment target -- must be a variable, index, or property expression", ctx);
    }

    private boolean isValidAssignmentTarget(Expression target) {
        // e.g. p = 123
        if( target instanceof VariableExpression && !isInsideParentheses(target) )
            return true;
        // e.g. obj.p = 123
        if( target instanceof PropertyExpression )
            return true;
        // e.g. map[a] = 123 OR map['a'] = 123 OR map["$a"] = 123
        if( target instanceof BinaryExpression be && be.getOperation().getType() == Types.LEFT_SQUARE_BRACKET )
            return true;
        return false;
    }

    private Statement expressionStatement(ExpressionStatementContext ctx) {
        var base = expression(ctx.expression());
        var expression = ctx.argumentList() != null
            ? methodCall(base, argumentList(ctx.argumentList()))
            : base;
        return ast( stmt(expression), ctx );
    }

    /// EXPRESSIONS

    private Expression expression(ExpressionContext ctx) {
        if( ctx instanceof AddSubExprAltContext asac )
            return ast( binary(asac.left, asac.op, asac.right), asac );

        if( ctx instanceof BitwiseAndExprAltContext baac )
            return ast( binary(baac.left, baac.op, baac.right), baac );

        if( ctx instanceof BitwiseOrExprAltContext boac )
            return ast( binary(boac.left, boac.op, boac.right), boac );

        if( ctx instanceof ConditionalExprAltContext cac )
            return ast( ternary(cac), cac );

        if( ctx instanceof EqualityExprAltContext eac )
            return ast( binary(eac.left, eac.op, eac.right), eac );

        if( ctx instanceof ExclusiveOrExprAltContext eoac )
            return ast( binary(eoac.left, eoac.op, eoac.right), eoac );

        if( ctx instanceof LogicalAndExprAltContext laac )
            return ast( binary(laac.left, laac.op, laac.right), laac );

        if( ctx instanceof LogicalOrExprAltContext loac )
            return ast( binary(loac.left, loac.op, loac.right), loac );

        if( ctx instanceof MultDivExprAltContext mdac )
            return ast( binary(mdac.left, mdac.op, mdac.right), mdac );

        if( ctx instanceof PathExprAltContext pac )
            return ast( pathExpression(pac), pac );

        if( ctx instanceof PowerExprAltContext pac )
            return ast( binary(pac.left, pac.op, pac.right), pac );

        if( ctx instanceof RegexExprAltContext rac )
            return ast( binary(rac.left, rac.op, rac.right), rac );

        if( ctx instanceof RelationalExprAltContext rac )
            return ast( binary(rac.left, rac.op, rac.right), rac );

        if( ctx instanceof RelationalCastExprAltContext rcac ) {
            var operand = expression(rcac.expression());
            var type = type(rcac.type());
            return ast( asX(type, operand), rcac );
        }

        if( ctx instanceof RelationalTypeExprAltContext rtac ) {
            var right = ast( new ClassExpression(type(rtac.type(), false)), rtac.type() );
            return ast( binary(rtac.left, rtac.op, right), rtac );
        }

        if( ctx instanceof ShiftExprAltContext sac )
            return ast( shift(sac), sac );

        if( ctx instanceof UnaryAddExprAltContext uaac )
            return ast( unaryAdd(expression(uaac.expression()), uaac.op), uaac );

        if( ctx instanceof UnaryNotExprAltContext unac )
            return ast( unaryNot(expression(unac.expression()), unac.op), unac );

        if( ctx instanceof IncompleteExprAltContext iac ) {
            var object = expression(iac.expression());
            var result = ast( propX(object, ""), iac );
            collectSyntaxError(new SyntaxException("Incomplete expression", result));
            return result;
        }

        throw createParsingFailedException("Invalid expression: " + ctx.getText(), ctx);
    }

    private Expression binary(ExpressionContext left, Token op, ExpressionContext right) {
        return binX(expression(left), token(op), expression(right));
    }

    private Expression binary(ExpressionContext left, Token op, Expression right) {
        return binX(expression(left), token(op), right);
    }

    private Expression shift(ShiftExprAltContext ctx) {
        var left = expression(ctx.left);
        var right = expression(ctx.right);

        if( ctx.riOp != null )
            return new RangeExpression(left, right, true);
        if( ctx.reOp != null )
            return new RangeExpression(left, right, false, true);

        org.codehaus.groovy.syntax.Token op = null;
        if( ctx.dlOp != null )
            op = token(ctx.dlOp, 2);
        if( ctx.dgOp != null )
            op = token(ctx.dgOp, 2);
        if( ctx.tgOp != null )
            op = token(ctx.tgOp, 3);

        return binX(left, op, right);
    }

    private Expression ternary(ConditionalExprAltContext ctx) {
        if( ctx.ELVIS() != null )
            return elvisX(expression(ctx.condition), expression(ctx.fb));

        var condition = ast( boolX(expression(ctx.condition)), ctx.condition );
        return ternaryX(condition, expression(ctx.tb), expression(ctx.fb));
    }

    private Expression unaryAdd(Expression expression, Token op) {
        if( op.getType() == ScriptParser.ADD )
            return new UnaryPlusExpression(expression);

        if( op.getType() == ScriptParser.SUB )
            return new UnaryMinusExpression(expression);

        throw new IllegalStateException();
    }

    private Expression unaryNot(Expression expression, Token op) {
        if( op.getType() == ScriptParser.NOT )
            return new NotExpression(expression);

        if( op.getType() == ScriptParser.BITNOT )
            return new BitwiseNegationExpression(expression);

        throw new IllegalStateException();
    }

    /// -- PATH EXPRESSIONS

    private Expression pathExpression(PathExprAltContext ctx) {
        try {
            var primary = primary(ctx.primary());
            return (Expression) ctx.pathElement().stream()
                .map(el -> (Object) el)
                .reduce(primary, (acc, el) -> pathElement((Expression) acc, (PathElementContext) el));
        }
        catch( IllegalStateException e ) {
            throw createParsingFailedException("Invalid expression: " + ctx.getText(), ctx);
        }
    }

    private Expression pathElement(Expression expression, PathElementContext ctx) {
        if( ctx instanceof PropertyPathExprAltContext pac )
            return ast( pathPropertyElement(expression, pac), expression, pac );

        if( ctx instanceof ClosurePathExprAltContext cac ) {
            var closure = closure(cac.closure());
            return ast( pathClosureElement(expression, closure), expression, cac );
        }

        if( ctx instanceof ClosureWithLabelsPathExprAltContext clac ) {
            var closure = closureWithLabels(clac.closureWithLabels());
            return ast( pathClosureElement(expression, closure), expression, clac );
        }

        if( ctx instanceof ArgumentsPathExprAltContext aac )
            return ast( pathArgumentsElement(expression, aac.arguments()), expression, aac );

        if( ctx instanceof IndexPathExprAltContext iac )
            return ast( pathIndexElement(expression, iac.indexPropertyArgs()), expression, iac );

        throw new IllegalStateException();
    }

    private Expression pathPropertyElement(Expression expression, PropertyPathExprAltContext ctx) {
        var property = namedProperty(ctx.namedProperty());
        var safe = ctx.SAFE_DOT() != null || ctx.SPREAD_DOT() != null;
        var result = new PropertyExpression(expression, property, safe);
        if( ctx.SPREAD_DOT() != null )
            result.setSpreadSafe(true);
        return result;
    }

    private Expression namedProperty(NamedPropertyContext ctx) {
        if( ctx.keywords() != null )
            return ast( constX(keywords(ctx.keywords())), ctx );

        if( ctx.identifier() != null )
            return ast( constX(identifier(ctx.identifier())), ctx );

        if( ctx.stringLiteral() != null )
            return ast( string(ctx.stringLiteral()), ctx );

        throw new IllegalStateException();
    }

    private Expression pathClosureElement(Expression expression, Expression closure) {
        if( expression instanceof MethodCallExpression mce ) {
            // normal arguments, e.g. 1, 2
            if ( !(mce.getArguments() instanceof ArgumentListExpression) )
                throw new IllegalStateException();

            var arguments = (ArgumentListExpression) mce.getArguments();
            arguments.addExpression(closure);
            return mce;
        }

        var arguments = ast( args(closure), closure );

        // e.g. obj.m { }
        if( expression instanceof PropertyExpression pe )
            return propMethodCall(pe, arguments);

        // e.g. m { }, "m" { }
        if( isConstMethodName(expression) )
            return thisMethodCall(expression, arguments);

        // e.g. <expr> { } -> <expr>.call { }
        return callMethodCall(expression, arguments);
    }

    private Expression pathArgumentsElement(Expression caller, ArgumentsContext ctx) {
        var arguments = argumentList(ctx.argumentList());
        if( ctx.COMMA() != null )
            arguments.putNodeMetaData(ASTNodeMarker.TRAILING_COMMA, Boolean.TRUE);
        return ast( methodCall(caller, arguments), caller, ctx );
    }

    private Expression pathIndexElement(Expression expression, IndexPropertyArgsContext ctx) {
        var elements = expressionList(ctx.expressionList());

        Expression index;
        if( elements.size() > 1 ) {
            // e.g. a[1, 2]
            var list = listX(elements);
            list.setWrapped(true);
            index = list;
        }
        else {
            // e.g. a[1]
            index = elements.get(0);
        }

        return indexX(expression, ast(index, ctx));
    }

    /// -- PRIMARY EXPRESSIONS

    private Expression primary(PrimaryContext ctx) {
        if( ctx instanceof IdentifierPrmrAltContext iac )
            return ast( variableName(iac.identifier()), iac );

        if( ctx instanceof LiteralPrmrAltContext lac )
            return ast( literal(lac.literal()), lac );

        if( ctx instanceof GstringPrmrAltContext gac )
            return ast( gstring(gac.gstring()), gac );

        if( ctx instanceof NewPrmrAltContext nac )
            return ast( creator(nac.creator()), nac );

        if( ctx instanceof ParenPrmrAltContext pac )
            return ast( parExpression(pac.parExpression()), pac );

        if( ctx instanceof ClosurePrmrAltContext cac )
            return ast( closure(cac.closure()), cac );

        if( ctx instanceof ListPrmrAltContext lac )
            return ast( list(lac.list()), lac );

        if( ctx instanceof MapPrmrAltContext mac )
            return ast( map(mac.map()), mac );

        if( ctx instanceof BuiltInTypePrmrAltContext bac )
            return ast( builtInType(bac.builtInType()), bac );

        throw createParsingFailedException("Invalid expression: " + ctx.getText(), ctx);
    }

    private Expression builtInType(BuiltInTypeContext ctx) {
        return varX(ctx.getText());
    }

    private String identifier(IdentifierContext ctx) {
        return ctx.getText();
    }

    private String keywords(KeywordsContext ctx) {
        return ctx.getText();
    }

    private Expression literal(LiteralContext ctx) {
        if( ctx instanceof IntegerLiteralAltContext iac )
            return ast( integerLiteral(iac), iac );

        if( ctx instanceof FloatingPointLiteralAltContext fac )
            return ast( floatingPointLiteral(fac), fac );

        if( ctx instanceof StringLiteralAltContext sac )
            return ast( string(sac.stringLiteral()), sac );

        if( ctx instanceof BooleanLiteralAltContext bac )
            return ast( constX("true".equals(bac.getText())), bac );

        if( ctx instanceof NullLiteralAltContext nac )
            return ast( constX(null), nac );

        throw createParsingFailedException("Invalid expression: " + ctx.getText(), ctx);
    }

    private Expression integerLiteral(IntegerLiteralAltContext ctx) {
        Number num = null;
        try {
            num = Numbers.parseInteger(ctx.getText());
        }
        catch( Exception e ) {
            numberFormatError = new Tuple2(ctx, e);
        }

        return constX(num, true);
    }

    private Expression floatingPointLiteral(FloatingPointLiteralAltContext ctx) {
        Number num = null;
        try {
            num = Numbers.parseDecimal(ctx.getText());
        }
        catch( Exception e ) {
            numberFormatError = new Tuple2(ctx, e);
        }

        return constX(num, true);
    }

    private ConstantExpression string(ParserRuleContext ctx) {
        var text = ctx.getText();
        var result = constX(stringLiteral(text));
        result.putNodeMetaData(ASTNodeMarker.VERBATIM_TEXT, text);
        return result;
    }

    private String stringLiteral(StringLiteralContext ctx) {
        return stringLiteral(ctx.getText());
    }

    private String stringLiteral(String text) {
        var startsWithSlash = text.startsWith(SLASH_STR);

        if( text.startsWith(TSQ_STR) || text.startsWith(TDQ_STR) ) {
            text = StringUtils.removeCR(text);
            text = StringUtils.trimQuotations(text, 3);
        }
        else if( text.startsWith(SQ_STR) || text.startsWith(DQ_STR) || startsWithSlash ) {
            // the slashy string can span rows, so we have to remove CR for it
            if( startsWithSlash )
                text = StringUtils.removeCR(text);
            text = StringUtils.trimQuotations(text, 1);
        }

        var slashyType = startsWithSlash
            ? StringUtils.SLASHY
            : StringUtils.NONE_SLASHY;

        return StringUtils.replaceEscapes(text, slashyType);
    }

    private Expression gstring(GstringContext ctx) {
        var text = ctx.getText();
        var beginQuotation = beginQuotation(text);
        var verbatimText = stringLiteral(text);
        var builder = new GStringBuilder();

        for( var part : ctx.gstringDqPart() ) {
            if( part instanceof GstringDqTextAltContext tac )
                builder.appendString(ast( gstringText(tac, beginQuotation), tac ));

            if( part instanceof GstringDqPathAltContext pac )
                builder.appendValue(ast( gstringPath(pac), pac ));

            if( part instanceof GstringDqExprAltContext eac )
                builder.appendValue(expression(eac.expression()));
        }

        for( var part : ctx.gstringTdqPart() ) {
            if( part instanceof GstringTdqTextAltContext tac )
                builder.appendString(ast( gstringText(tac, beginQuotation), tac ));

            if( part instanceof GstringTdqPathAltContext pac )
                builder.appendValue(ast( gstringPath(pac), pac ));

            if( part instanceof GstringTdqExprAltContext eac )
                builder.appendValue(expression(eac.expression()));
        }

        var result = builder.build(verbatimText);
        result.putNodeMetaData(ASTNodeMarker.QUOTE_CHAR, beginQuotation);
        return result;
    }

    /**
     * Builder for GStringExpression that inserts empty strings
     * to ensure that there are n+1 strings for n values.
     * 
     * @see org.codehaus.groovy.runtime.GStringUtil.writeToImpl()
     */
    private static class GStringBuilder {
        private final List<ConstantExpression> strings = new ArrayList<>();
        private final List<Expression> values = new ArrayList<>();
        private boolean appendEmptyString = true;

        public void appendString(ConstantExpression string) {
            strings.add(string);
            appendEmptyString = false;
        }

        public void appendValue(Expression value) {
            if( appendEmptyString )
                appendString(constX(""));
            values.add(value);
            appendEmptyString = true;
        }

        public GStringExpression build(String verbatimText) {
            return new GStringExpression(verbatimText, strings, values);
        }
    }

    private String beginQuotation(String text) {
        if( text.startsWith(TDQ_STR) )
            return TDQ_STR;
        if( text.startsWith(DQ_STR) )
            return DQ_STR;

        throw new IllegalStateException();
    }

    private ConstantExpression gstringText(ParserRuleContext ctx, String beginQuotation) {
        var text = ctx.getText();
        var quotedText = new StringBuilder(text)
            .insert(0, beginQuotation)
            .append(beginQuotation)
            .toString();
        var result = constX(stringLiteral(quotedText));
        result.putNodeMetaData(ASTNodeMarker.VERBATIM_TEXT, text);
        return result;
    }

    private Expression gstringPath(ParserRuleContext ctx) {
        var names = ctx.getText().split("\\.");
        int currentLine = ctx.getStart().getLine();
        int currentChar = ctx.getStart().getCharPositionInLine() + 1;
        var varName = names[0].substring(1);
        Expression result = varX(varName);
        currentChar += 1;
        result.setLineNumber(currentLine);
        result.setColumnNumber(currentChar);
        currentChar += varName.length();
        result.setLastLineNumber(currentLine);
        result.setLastColumnNumber(currentChar);

        for( int i = 1; i < names.length; i++ ) {
            var propName = names[i];
            var property = constX(propName);
            currentChar += 1;
            property.setLineNumber(currentLine);
            property.setColumnNumber(currentChar);
            currentChar += propName.length();
            property.setLastLineNumber(currentLine);
            property.setLastColumnNumber(currentChar);
            result = ast( propX(result, property), result, property );
        }

        return result;
    }

    private Expression creator(CreatorContext ctx) {
        var type = createdName(ctx.createdName());
        var arguments = argumentList(ctx.arguments().argumentList());
        if( ctx.arguments().COMMA() != null )
            arguments.putNodeMetaData(ASTNodeMarker.TRAILING_COMMA, Boolean.TRUE);
        return ctorX(type, arguments);
    }

    private Expression parExpression(ParExpressionContext ctx) {
        var expression = expression(ctx.expression());
        expression.getNodeMetaData(ASTNodeMarker.INSIDE_PARENTHESES_LEVEL, k -> new AtomicInteger()).getAndAdd(1);
        return expression;
    }

    private Expression closure(ClosureContext ctx) {
        var params = formalParameterList(ctx.formalParameterList());
        var code = blockStatements(ctx.blockStatements());
        return ast( closureX(params, code), ctx );
    }

    private Expression closureWithLabels(ClosureWithLabelsContext ctx) {
        var params = formalParameterList(ctx.formalParameterList());
        var code = blockStatementsWithLabels(ctx.blockStatementsWithLabels());
        return ast( closureX(params, code), ctx );
    }

    private BlockStatement blockStatementsWithLabels(BlockStatementsWithLabelsContext ctx) {
        var statements = ctx.statementOrLabeled().stream()
            .map(this::statementOrLabeled)
            .toList();
        return ast( block(new VariableScope(), statements), ctx );
    }

    private Statement statementOrLabeled(StatementOrLabeledContext ctx) {
        if( ctx.identifier() != null ) {
            var label = identifier(ctx.identifier());
            var result = statementOrLabeled(ctx.statementOrLabeled());
            result.addStatementLabel(label);
            return result;
        }
        else {
            return statement(ctx.statement());
        }
    }

    private Expression list(ListContext ctx) {
        if( ctx.COMMA() != null && ctx.expressionList() == null )
            throw createParsingFailedException("Empty list literal should not contain any comma(,)", ctx.COMMA());

        var result = listX(expressionList(ctx.expressionList()));
        if( ctx.COMMA() != null )
            result.putNodeMetaData(ASTNodeMarker.TRAILING_COMMA, Boolean.TRUE);
        return result;
    }

    private List<Expression> expressionList(ExpressionListContext ctx) {
        if( ctx == null )
            return Collections.emptyList();
        
        return ctx.expression().stream()
            .map(this::expression)
            .toList();
    }

    private Expression map(MapContext ctx) {
        if( ctx.mapEntryList() == null )
            return new MapExpression();

        var entries = ctx.mapEntryList().mapEntry().stream()
            .map(this::mapEntry)
            .toList();
        var result = mapX(entries);
        if( ctx.COMMA() != null )
            result.putNodeMetaData(ASTNodeMarker.TRAILING_COMMA, Boolean.TRUE);
        return result;
    }

    private MapEntryExpression mapEntry(MapEntryContext ctx) {
        var key = mapEntryLabel(ctx.mapEntryLabel());
        var value = expression(ctx.expression());
        return ast( entryX(key, value), ctx );
    }

    private Expression mapEntryLabel(MapEntryLabelContext ctx) {
        if( ctx.keywords() != null )
            return ast( constX(keywords(ctx.keywords())), ctx );

        if( ctx.primary() != null ) {
            var expression = primary(ctx.primary());
            return expression instanceof VariableExpression ve && !isInsideParentheses(ve)
                ? ast( constX(ve.getName()), ve )
                : ast( expression, ctx );
        }

        throw createParsingFailedException("Unsupported map entry label: " + ctx.getText(), ctx);
    }

    private Expression methodCall(Expression caller, Expression arguments) {
        // e.g. (obj.x)(), (obj.@x)()
        if( isInsideParentheses(caller) )
            return callMethodCall(caller, arguments);

        // e.g. obj.a(1, 2)
        if( caller instanceof PropertyExpression pe )
            return propMethodCall(pe, arguments);

        // e.g. m(), "m"()
        if( isConstMethodName(caller) )
            return thisMethodCall(caller, arguments);

        // e.g. <expr>(<args>) -> <expr>.call(<args>)
        return callMethodCall(caller, arguments);
    }

    private boolean isConstMethodName(Expression caller) {
        if( caller instanceof VariableExpression )
            return true;
        if( caller instanceof ConstantExpression ce && ce.getValue() instanceof String )
            return true;
        return false;
    }

    private Expression propMethodCall(PropertyExpression caller, Expression arguments) {
        var result = callX(caller.getObjectExpression(), caller.getProperty(), arguments);
        result.setImplicitThis(false);
        result.setSafe(caller.isSafe());
        result.setSpreadSafe(caller.isSpreadSafe());

        // method call obj*.m() -> safe=false and spreadSafe=true
        // property access obj*.p -> safe=true and spreadSafe=true
        if( caller.isSpreadSafe() )
            result.setSafe(false);

        return ast( result, caller, arguments );
    }

    private Expression thisMethodCall(Expression caller, Expression arguments) {
        var object = THIS_EXPRESSION;
        object.setColumnNumber(caller.getColumnNumber());
        object.setLineNumber(caller.getLineNumber());

        var name = caller instanceof VariableExpression ve
            ? ast( constX(ve.getText()), ve )
            : caller;

        return ast( callX(object, name, arguments), caller, arguments );
    }

    private Expression callMethodCall(Expression caller, Expression arguments) {
        var call = callX(caller, CALL_STR, arguments);
        call.setImplicitThis(false);
        return ast( call, caller, arguments );
    }

    private Expression argumentList(ArgumentListContext ctx) {
        if( ctx == null )
            return new ArgumentListExpression();

        var arguments = new ArrayList<Expression>();
        var namedArgs = new ArrayList<MapEntryExpression>();

        for( var el : ctx.argumentListElement() ) {
            if( el.expression() != null ) {
                arguments.add(expression(el.expression()));
            }
            else if( el.namedArg() != null ) {
                var namedArg = namedArg(el.namedArg());
                checkDuplicateNamedArg(namedArgs, namedArg);
                namedArgs.add(namedArg);
            }
        }

        if( !namedArgs.isEmpty() )
            arguments.add(0, new NamedArgumentListExpression(namedArgs));

        return ast( args(arguments), ctx );
    }

    private MapEntryExpression namedArg(NamedArgContext ctx) {
        var key = namedProperty(ctx.namedProperty());
        var value = expression(ctx.expression());
        return ast( new MapEntryExpression(key, value), ctx );
    }

    private void checkDuplicateNamedArg(List<MapEntryExpression> namedArgs, MapEntryExpression namedArg) {
        var name = namedArg.getKeyExpression().getText();
        for( var arg : namedArgs ) {
            if( arg.getKeyExpression().getText().equals(name) )
                throw createParsingFailedException("Duplicated named argument '" + name + "' found", namedArg);
        }
    }

    /// MISCELLANEOUS

    private Parameter[] formalParameterList(FormalParameterListContext ctx) {
        // NOTE: implicit `it` parameter is deprecated, but allow it for now
        if( ctx == null )
            return Parameter.EMPTY_ARRAY;

        var params = ctx.formalParameter().stream()
            .map(this::formalParameter)
            .toList();
        for( int n = params.size(), i = n - 1; i >= 0; i -= 1 ) {
            var param = params.get(i);
            for( var other : params ) {
                if( other == param )
                    continue;
                if( other.getName().equals(param.getName()) )
                    throw createParsingFailedException("Duplicated parameter '" + param.getName() + "' found", param);
            }
        }

        return params.toArray(Parameter.EMPTY_ARRAY);
    }

    private Parameter formalParameter(FormalParameterContext ctx) {
        var type = legacyType(ctx.legacyType());
        var name = identifier(ctx.identifier());
        var defaultValue = ctx.expression() != null
            ? expression(ctx.expression())
            : null;
        var result = ast( param(type, name, defaultValue), ctx );
        checkInvalidVarName(name, result);
        result.putNodeMetaData("_START_NAME", tokenPosition(ctx.identifier()));
        return result;
    }

    private org.codehaus.groovy.syntax.Token token(Token token) {
        return token(token, 1);
    }

    private org.codehaus.groovy.syntax.Token token(Token token, int cardinality) {
        var tokenText = token.getText();
        var tokenType = token.getType();
        var text = cardinality == 1
            ? tokenText
            : tokenText.repeat(cardinality);
        var type = tokenType == RANGE_EXCLUSIVE_RIGHT || tokenType == RANGE_INCLUSIVE
            ? Types.RANGE_OPERATOR
            : Types.lookup(text, Types.ANY);
        return new org.codehaus.groovy.syntax.Token( type, text, token.getLine(), token.getCharPositionInLine() + 1 );
    }

    private ClassNode createdName(CreatedNameContext ctx) {
        if( ctx.qualifiedClassName() != null ) {
            var classNode = qualifiedClassName(ctx.qualifiedClassName());
            if( ctx.typeArguments() != null )
                classNode.setGenericsTypes( typeArguments(ctx.typeArguments()) );
            return classNode;
        }

        if( ctx.primitiveType() != null )
            return primitiveType(ctx.primitiveType());

        throw createParsingFailedException("Unrecognized created name: " + ctx.getText(), ctx);
    }

    private ClassNode primitiveType(PrimitiveTypeContext ctx) {
        var classNode = ClassHelper.make(ctx.getText()).getPlainNodeReference(false);
        return ast( classNode, ctx );
    }

    private ClassNode qualifiedClassName(QualifiedClassNameContext ctx) {
        return qualifiedClassName(ctx, true);
    }

    private ClassNode qualifiedClassName(QualifiedClassNameContext ctx, boolean allowProxy) {
        var text = ctx.getText();
        var classNode = ClassHelper.make(text);
        if( text.contains(".") )
            classNode.putNodeMetaData(ASTNodeMarker.FULLY_QUALIFIED, true);

        if( classNode.isUsingGenerics() && allowProxy ) {
            var proxy = ClassHelper.makeWithoutCaching(classNode.getName());
            proxy.setRedirect(classNode);
            return proxy;
        }

        return ast( classNode, ctx );
    }

    private ClassNode type(TypeContext ctx) {
        return type(ctx, true);
    }

    private ClassNode type(TypeContext ctx, boolean allowProxy) {
        if( ctx == null )
            return ClassHelper.dynamicType();

        if( ctx.qualifiedClassName() != null ) {
            var classNode = qualifiedClassName(ctx.qualifiedClassName(), allowProxy);
            if( ctx.typeArguments() != null )
                classNode.setGenericsTypes( typeArguments(ctx.typeArguments()) );
            return classNode;
        }

        if( ctx.primitiveType() != null )
            return primitiveType(ctx.primitiveType());

        throw createParsingFailedException("Unrecognized type: " + ctx.getText(), ctx);
    }

    private GenericsType[] typeArguments(TypeArgumentsContext ctx) {
        return ctx.type().stream()
            .map(this::genericsType)
            .toArray(GenericsType[]::new);
    }

    private GenericsType genericsType(TypeContext ctx) {
        return ast( new GenericsType(type(ctx)), ctx );
    }

    private ClassNode legacyType(ParserRuleContext ctx) {
        var result = ClassHelper.dynamicType();
        if( ctx != null )
            result.putNodeMetaData(ASTNodeMarker.LEGACY_TYPE, ctx.getText());
        return result;
    }

    /// COMMENTS

    private void saveLeadingComments(ASTNode node, ParserRuleContext ctx) {
        var comments = new ArrayList<String>();
        var child = ctx;
        while( saveLeadingComments0(child, comments) )
            child = child.getParent();

        if( !comments.isEmpty() )
            node.putNodeMetaData(ASTNodeMarker.LEADING_COMMENTS, comments);
    }

    private boolean saveLeadingComments0(ParserRuleContext ctx, List<String> comments) {
        var parent = ctx.getParent();
        if( parent == null )
            return false;

        // find index of token among siblings
        var siblings = parent.children;
        int i = 0;
        while( i < siblings.size() && siblings.get(i) != ctx ) {
            i++;
        }

        // check parent context for additional comments
        if( i == 0 )
            return true;

        // prepend each comment/newline to node
        var added = false;
        for( int j = i - 1; j >= 0; j-- ) {
            var sibling = siblings.get(j);
            if( !(sibling instanceof NlsContext || sibling instanceof SepContext) )
                break;

            var newlines = sibling instanceof NlsContext
                ? ((NlsContext) sibling).NL()
                : ((SepContext) sibling).NL();

            for( int k = newlines.size() - 1; k >= 0; k-- ) {
                var text = newlines.get(k).getText();
                comments.add(text);
                added = true;
            }
        }

        // remove leading newline
        if( added ) {
            var last = comments.get(comments.size() - 1);
            if( "\n".equals(last) ) {
                comments.remove(comments.size() - 1);
            }
            else if( last.startsWith("#!") ) {
                moduleNode.setShebang(last);
                comments.remove(comments.size() - 1);
            }
        }

        return false;
    }

    private void saveTrailingComment(ASTNode node, ParserRuleContext ctx) {
        var child = ctx;
        while( saveTrailingComment0(child, node) )
            child = child.getParent();
    }

    private boolean saveTrailingComment0(ParserRuleContext ctx, ASTNode node) {
        var parent = ctx.getParent();
        if( parent == null )
            return false;

        // find index of token among siblings
        var siblings = parent.children;
        int i = 0;
        while( i < siblings.size() && siblings.get(i) != ctx ) {
            i++;
        }

        // check parent context for additional comments
        if( i == siblings.size() - 1 )
            return true;

        // save trailing inline comment
        var next = siblings.get(i + 1);
        if( next instanceof SepContext sep ) {
            var text = sep.children.get(0).getText();
            if( !";".equals(text) && !"\n".equals(text) )
                node.putNodeMetaData(ASTNodeMarker.TRAILING_COMMENT, text);
        }

        return false;
    }

    /// HELPERS

    private boolean isInsideParentheses(NodeMetaDataHandler nodeMetaDataHandler) {
        Number insideParenLevel = nodeMetaDataHandler.getNodeMetaData(ASTNodeMarker.INSIDE_PARENTHESES_LEVEL);
        return insideParenLevel != null && insideParenLevel.intValue() > 0;
    }

    private CompilationFailedException createParsingFailedException(String msg, ParserRuleContext ctx) {
        return createParsingFailedException(
            new SyntaxException(msg,
                ctx.start.getLine(),
                ctx.start.getCharPositionInLine() + 1,
                ctx.stop.getLine(),
                ctx.stop.getCharPositionInLine() + 1 + ctx.stop.getText().length()));
    }

    private CompilationFailedException createParsingFailedException(String msg, Tuple2<Integer,Integer> start, Tuple2<Integer,Integer> end) {
        return createParsingFailedException(
            new SyntaxException(msg,
                start.getV1(),
                start.getV2(),
                end.getV1(),
                end.getV2()));
    }

    private CompilationFailedException createParsingFailedException(String msg, ASTNode node) {
        return createParsingFailedException(
            new SyntaxException(msg,
                node.getLineNumber(),
                node.getColumnNumber(),
                node.getLastLineNumber(),
                node.getLastColumnNumber()));
    }

    private CompilationFailedException createParsingFailedException(String msg, TerminalNode node) {
        return createParsingFailedException(msg, node.getSymbol());
    }

    private CompilationFailedException createParsingFailedException(String msg, Token token) {
        return createParsingFailedException(
            new SyntaxException(msg,
                token.getLine(),
                token.getCharPositionInLine() + 1,
                token.getLine(),
                token.getCharPositionInLine() + 1 + token.getText().length()));
    }

    private CompilationFailedException createParsingFailedException(Throwable t) {
        if( t instanceof SyntaxException se )
            collectSyntaxError(se);

        else if( t instanceof GroovySyntaxError gse )
            collectSyntaxError(
                    new SyntaxException(
                            gse.getMessage(),
                            gse,
                            gse.getLine(),
                            gse.getColumn()));

        else if( t instanceof Exception e )
            collectException(e);

        return new CompilationFailedException(
                CompilePhase.PARSING.getPhaseNumber(),
                sourceUnit,
                t);
    }

    private void collectSyntaxError(SyntaxException e) {
        sourceUnit.getErrorCollector().addFatalError(new SyntaxErrorMessage(e, sourceUnit));
    }

    private void collectException(Exception e) {
        sourceUnit.getErrorCollector().addException(e, sourceUnit);
    }

    private void collectWarning(String text, ASTNode node) {
        sourceUnit.addWarning(text, node);
    }

    private void removeErrorListeners() {
        lexer.removeErrorListeners();
        parser.removeErrorListeners();
    }

    private void addErrorListeners() {
        lexer.addErrorListener(createANTLRErrorListener());
        parser.addErrorListener(createANTLRErrorListener());
    }

    private ANTLRErrorListener createANTLRErrorListener() {
        return new ANTLRErrorListener() {
            @Override
            public void syntaxError(Recognizer recognizer, Object offendingSymbol, int line, int charPositionInLine, String msg, RecognitionException e) {
                collectSyntaxError(new SyntaxException(msg, line, charPositionInLine + 1));
            }

            @Override
            public void reportAmbiguity(Parser recognizer, DFA dfa, int startIndex, int stopIndex, boolean exact, BitSet ambigAlts, ATNConfigSet configs) {}

            @Override
            public void reportAttemptingFullContext(Parser recognizer, DFA dfa, int startIndex, int stopIndex, BitSet conflictingAlts, ATNConfigSet configs) {}

            @Override
            public void reportContextSensitivity(Parser recognizer, DFA dfa, int startIndex, int stopIndex, int prediction, ATNConfigSet configs) {}
        };
    }

    private static final String CALL_STR = "call";
    private static final String SLASH_STR = "/";
    private static final String TDQ_STR = "\"\"\"";
    private static final String TSQ_STR = "'''";
    private static final String SQ_STR = "'";
    private static final String DQ_STR = "\"";

    private static final List<String> GROOVY_KEYWORDS = List.of(
        Types.getText(Types.KEYWORD_ABSTRACT),
        Types.getText(Types.KEYWORD_ASSERT),
        Types.getText(Types.KEYWORD_BREAK),
        Types.getText(Types.KEYWORD_CASE),
        Types.getText(Types.KEYWORD_CLASS),
        Types.getText(Types.KEYWORD_CONST),
        Types.getText(Types.KEYWORD_CONTINUE),
        Types.getText(Types.KEYWORD_DEFAULT),
        Types.getText(Types.KEYWORD_DO),
        Types.getText(Types.KEYWORD_EXTENDS),
        Types.getText(Types.KEYWORD_FINAL),
        Types.getText(Types.KEYWORD_FINALLY),
        Types.getText(Types.KEYWORD_FOR),
        Types.getText(Types.KEYWORD_GOTO),
        Types.getText(Types.KEYWORD_IMPLEMENTS),
        // Types.getText(Types.KEYWORD_IMPORT),
        Types.getText(Types.KEYWORD_INTERFACE),
        Types.getText(Types.KEYWORD_NATIVE),
        Types.getText(Types.KEYWORD_PACKAGE),
        Types.getText(Types.KEYWORD_PRIVATE),
        Types.getText(Types.KEYWORD_PROTECTED),
        Types.getText(Types.KEYWORD_PUBLIC),
        Types.getText(Types.KEYWORD_STATIC),
        Types.getText(Types.KEYWORD_SUPER),
        Types.getText(Types.KEYWORD_SWITCH),
        Types.getText(Types.KEYWORD_SYNCHRONIZED),
        Types.getText(Types.KEYWORD_THIS),
        Types.getText(Types.KEYWORD_THROWS),
        Types.getText(Types.KEYWORD_TRANSIENT),
        Types.getText(Types.KEYWORD_VOID),
        Types.getText(Types.KEYWORD_VOLATILE),
        Types.getText(Types.KEYWORD_WHILE)
    );

}
