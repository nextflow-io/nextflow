/*
 * Copyright 2024-2026, Seqera Labs
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
package nextflow.script.control;

import java.lang.reflect.Modifier;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import nextflow.script.ast.AgentNode;
import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.ast.AssignmentExpression;
import nextflow.script.ast.FeatureFlagNode;
import nextflow.script.ast.FunctionNode;
import nextflow.script.ast.OutputNode;
import nextflow.script.ast.ProcessNode;
import nextflow.script.ast.ProcessNodeV2;
import nextflow.script.ast.RecordNode;
import nextflow.script.ast.ScriptNode;
import nextflow.script.ast.ScriptVisitorSupport;
import nextflow.script.ast.WorkflowNode;
import nextflow.script.dsl.Namespace;
import nextflow.script.dsl.Ops;
import nextflow.script.types.Channel;
import nextflow.script.types.ParamsMap;
import nextflow.script.types.Record;
import nextflow.script.types.Tuple;
import nextflow.script.dsl.Types;
import nextflow.script.types.Value;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.FieldNode;
import org.codehaus.groovy.ast.GenericsType;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.PropertyNode;
import org.codehaus.groovy.ast.Variable;
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
import org.codehaus.groovy.ast.expr.ListExpression;
import org.codehaus.groovy.ast.expr.MapExpression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.PropertyExpression;
import org.codehaus.groovy.ast.expr.RangeExpression;
import org.codehaus.groovy.ast.expr.TernaryExpression;
import org.codehaus.groovy.ast.expr.TupleExpression;
import org.codehaus.groovy.ast.expr.UnaryMinusExpression;
import org.codehaus.groovy.ast.expr.UnaryPlusExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.ErrorCollector;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;
import org.codehaus.groovy.control.messages.WarningMessage;
import org.codehaus.groovy.syntax.Token;

import static nextflow.script.ast.ASTUtils.*;
import static nextflow.script.types.TypeCheckingUtils.*;
import static org.codehaus.groovy.syntax.Types.*;

/**
 * Resolve and validate the types of expressions.
 *
 * @see org.codehaus.groovy.transform.stc.StaticTypeCheckingVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class TypeCheckingVisitor extends ScriptVisitorSupport {

    private static final ClassNode CHANNEL_TYPE = ClassHelper.makeCached(Channel.class);
    private static final ClassNode PATH_TYPE = ClassHelper.makeCached(Path.class);
    private static final ClassNode PARAMS_TYPE = ClassHelper.makeCached(ParamsMap.class);
    private static final ClassNode RECORD_TYPE = ClassHelper.makeCached(Record.class);
    private static final ClassNode TUPLE_TYPE = ClassHelper.makeCached(Tuple.class);
    private static final ClassNode VALUE_TYPE = ClassHelper.makeCached(Value.class);

    private SourceUnit sourceUnit;

    private ErrorCollector errorCollector;

    public TypeCheckingVisitor(SourceUnit sourceUnit) {
        this(sourceUnit, sourceUnit.getErrorCollector());
    }

    /**
     * Create a type checker that reports to the given error collector, so that
     * type errors can be collected without failing a compilation.
     *
     * @param sourceUnit
     * @param errorCollector
     */
    public TypeCheckingVisitor(SourceUnit sourceUnit, ErrorCollector errorCollector) {
        this.sourceUnit = sourceUnit;
        this.errorCollector = errorCollector;
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    public void visit() {
        var moduleNode = sourceUnit.getAST();
        if( !(moduleNode instanceof ScriptNode sn) )
            return;
        for( var featureFlag : sn.getFeatureFlags() )
            visitFeatureFlag(featureFlag);
        if( sn.getParams() != null )
            visitParams(sn.getParams());
        for( var functionNode : sn.getFunctions() )
            visitFunction(functionNode);
        for( var processNode : sn.getProcesses() )
            visitProcess(processNode);
        for( var agentNode : sn.getAgents() )
            visitAgent(agentNode);
        for( var workflowNode : sn.getWorkflows() )
            visitWorkflow(workflowNode);
        if( sn.getOutputs() != null )
            visitOutputs(sn.getOutputs());
    }

    // script declarations

    @Override
    public void visitFeatureFlag(FeatureFlagNode node) {
        var fn = node.target;
        if( fn == null )
            return;
        var expectedType = fn.getType();
        var actualType = node.value.getType();
        if( !Types.isAssignableFrom(expectedType, actualType) )
            addError("Feature flag '" + node.name + "' expects a " + Types.getName(expectedType) + " but received a " + Types.getName(actualType), node);
    }

    private boolean hasParamsBlock;

    @Override
    public void visitParam(Parameter node) {
        visitParameter(node, true);
        hasParamsBlock = true;
    }

    private void visitParameter(Parameter node, boolean allowPathCoercion) {
        if( !node.hasInitialExpression() )
            return;
        var expectedType = node.getType();
        var actualType = getType(node.getInitialExpression());
        if( Types.isAssignableFrom(expectedType, actualType) )
            return;
        if( allowPathCoercion && Types.isEqual(expectedType, PATH_TYPE) && Types.isEqual(actualType, ClassHelper.STRING_TYPE) )
            return;
        addError("Parameter '" + node.getName() + "' with type " + Types.getName(expectedType) + " cannot be assigned to default value with type " + Types.getName(actualType), node);
    }

    private WorkflowNode currentWorkflow;

    @Override
    public void visitWorkflow(WorkflowNode node) {
        currentWorkflow = node;

        visit(node.main);

        if( node.isEntry() ) {
            visitWorkflowOutputs(node);
        }
        else {
            checkSingleNamedOutput(node.emits, "emit");
            visit(node.emits);
        }

        visit(node.onComplete);
        visit(node.onError);

        currentWorkflow = null;
    }

    private void visitWorkflowOutputs(WorkflowNode node) {
        var sn = (ScriptNode) sourceUnit.getAST();
        var outputs = sn.getOutputs();
        if( outputs != null && node.publishers.isEmpty() )
            addError("Output block is defined but entry workflow does not define a `publish:` section", node);
        if( outputs == null && !node.publishers.isEmpty() )
            addError("Entry workflow defines a `publish:` section but no output block is defined", node.publishers);
        if( outputs == null )
            return;

        for( var publisher : asBlockStatements(node.publishers) ) {
            var es = (ExpressionStatement)publisher;
            var ae = (AssignmentExpression)es.getExpression();

            var target = (VariableExpression)ae.getLeftExpression();
            var decl = outputs.declarations.stream()
                .filter(output -> output.getName().equals(target.getName()))
                .findFirst().orElse(null);
            if( decl == null ) {
                addError("Workflow output '" + target.getName() + "' was assigned in the entry workflow but not declared in the output block", publisher);
                continue;
            }
            target.setAccessedVariable(decl);

            var source = ae.getRightExpression();
            visit(source);
            var sourceType = getType(source);
            var targetType = asDataflowType(target.getType(), sourceType);
            if( !Types.isAssignableFrom(targetType, sourceType) )
                addError("Workflow output '" + target.getName() + "' with type " + Types.getName(targetType) + " cannot be assigned to value with type " + Types.getName(sourceType), ae);
        }
    }

    private static ClassNode asDataflowType(ClassNode type, ClassNode sourceType) {
        if( CHANNEL_TYPE.equals(type) || VALUE_TYPE.equals(type) )
            return type;
        if( VALUE_TYPE.equals(sourceType) )
            return makeType(VALUE_TYPE, type);
        return type;
    }

    @Override
    public void visitProcessV2(ProcessNodeV2 node) {
        visitProcessDirectives(node.directives);
        visit(node.stagers);
        if( !(node.when instanceof EmptyExpression) )
            addSoftError("Process `when` section is discouraged with static typing -- use conditional logic in the calling workflow instead", node.when);
        visit(node.when);
        visit(node.exec);
        visit(node.stub);
        checkSingleNamedOutput(node.outputs, "output");
        visit(node.outputs);
        visitProcessTopics(node.topics);
    }

    @Override
    public void visitAgent(AgentNode node) {
        visitProcessDirectives(node.directives);
        visit(node.outputs);
        visit(node.prompt);
    }

    private void visitProcessDirectives(Statement block) {
        asDirectives(block).forEach((call) -> {
            var arguments = asMethodCallArguments(call);
            if( arguments.size() == 1 && arguments.get(0) instanceof ClosureExpression ) {
                addWarning("Closure is no longer required for dynamic process directives", call.getMethodAsString(), call);
                return;
            }
            visit(call);
        });
    }

    private void checkSingleNamedOutput(Statement block, String section) {
        var outputs = asBlockStatements(block);
        if( outputs.size() != 1 )
            return;
        var output = ((ExpressionStatement) outputs.get(0)).getExpression();
        if( output instanceof AssignmentExpression ae )
            addError("Name should be omitted for a single " + section, ae);
    }

    private void visitProcessTopics(Statement block) {
        for( var stmt : asBlockStatements(block) ) {
            var es = (ExpressionStatement)stmt;
            var be = (BinaryExpression)es.getExpression();
            var source = be.getLeftExpression();
            var target = be.getRightExpression();
            visit(source);
            visit(target);
            var targetType = getType(target);
            if( !Types.isEqual(ClassHelper.STRING_TYPE, targetType) )
                addError("Topic name should be a String but was specified as a " + Types.getName(targetType), target);
        }
    }

    @Override
    public void visitFunction(FunctionNode node) {
        // visit parameters and code
        for( var parameter : node.getParameters() ) {
            checkLegacyType(parameter.getType(), parameter);
            visitParameter(parameter, false);
        }
        visit(node.getCode());

        // check return statements against declared return type
        var visitor = new ReturnStatementVisitor(sourceUnit, errorCollector);
        visitor.visit(node.getReturnType(), node.getCode());

        var inferredReturnType = visitor.getInferredReturnType();
        if( inferredReturnType != null && ClassHelper.isDynamicTyped(node.getReturnType()) )
            node.putNodeMetaData(ASTNodeMarker.INFERRED_RETURN_TYPE, inferredReturnType);
    }

    @Override
    public void visitOutput(OutputNode node) {
        var type = node.getType();
        var elementType = dataflowElementType(type);
        for( var stmt : asBlockStatements(node.body) ) {
            var call = asMethodCallX(stmt);
            if( checkPublishStatements(call, elementType) )
                continue;
            visitMethodCallExpression(call);
        }
    }

    private boolean checkPublishStatements(MethodCallExpression node, ClassNode elementType) {
        if( !"path".equals(node.getMethodAsString()) )
            return false;
        var args = asMethodCallArguments(node);
        if( args.size() != 1 )
            return false;
        var firstArg = args.get(0);
        if( !(firstArg instanceof ClosureExpression closure) )
            return false;
        if( closure.getParameters().length == 1 )
            closure.getParameters()[0].putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, elementType);
        var code = (BlockStatement) closure.getCode();
        boolean checkPublish = code.getStatements().stream().anyMatch(stmt -> (
            asPublishStatement(stmt) != null
        ));
        for( var stmt : code.getStatements() ) {
            if( checkPublish )
                checkPublishStatement(stmt);
            else
                visit(stmt);
        }
        return true;
    }

    private void checkPublishStatement(Statement node) {
        var publish = asPublishStatement(node);
        if( publish == null ) {
            addError("Statement is not a valid publish statement (hint: make sure right-hand side of `>>` expression is wrapped in parentheses)", node);
            return;
        }
        var source = publish.getLeftExpression();
        var target = publish.getRightExpression();
        visit(source);
        visit(target);
        var sourceType = getType(source);
        if( !isPathOrCollection(sourceType) )
            addError("Publish source should be a Path or Iterable<Path> but was specified as a " + Types.getName(sourceType), source);
        var targetType = getType(target);
        if( !Types.isAssignableFrom(ClassHelper.STRING_TYPE, targetType) )
            addError("Publish target should be a String but was specified as a " + Types.getName(targetType), target);
    }

    private BinaryExpression asPublishStatement(Statement node) {
        return node instanceof ExpressionStatement es
            && es.getExpression() instanceof BinaryExpression be
            && be.getOperation().getType() == RIGHT_SHIFT ? be : null;
    }

    private boolean isPathOrCollection(ClassNode type) {
        return Types.isAssignableFrom(PATH_TYPE, type)
            || Types.isAssignableFrom(ClassHelper.ITERABLE_TYPE, type);
    }

    // statements

    @Override
    public void visitExpressionStatement(ExpressionStatement node) {
        var exp = node.getExpression();
        if( exp instanceof AssignmentExpression ae ) {
            applyAssignment(ae);
            return;
        }
        super.visitExpressionStatement(node);
    }

    // expressions

    @Override
    public void visitDeclarationExpression(DeclarationExpression node) {
        checkLegacyType(node.getLeftExpression());
        applyAssignment(node);
    }

    private void checkLegacyType(Expression target) {
        if( target instanceof VariableExpression ve ) {
            checkLegacyType(ve.getType(), ve);
        }
        else if( target instanceof TupleExpression te ) {
            for( var el : te.getExpressions() )
                checkLegacyType(el);
        }
    }

    private void checkLegacyType(ClassNode type, ASTNode node) {
        if( type.getNodeMetaData(ASTNodeMarker.LEGACY_TYPE) != null )
            addError("Legacy type annotations are not supported with static typing -- use `name: Type` instead", node);
    }

    private void applyAssignment(BinaryExpression node) {
        var target = node.getLeftExpression();
        var source = node.getRightExpression();
        if( source instanceof EmptyExpression )
            return;
        visit(target);
        visit(source);
        var targetType = getType(target);
        var sourceType = getType(source);

        // a call with no return value (e.g. a function with a void return type or a
        // process/workflow that declares no outputs) cannot be assigned to a variable
        if( ClassHelper.isPrimitiveVoid(sourceType) ) {
            addError("Cannot assign the result of a call that does not return a value", node);
            return;
        }

        // a compound assignment (e.g. `total += 1`) is equivalent to applying the
        // corresponding binary operator and assigning the result (e.g. `total = total + 1`),
        // so resolve the operator result type and use it as the source type
        var opMethod = compoundAssignmentMethod(node.getOperation().getType());
        if( opMethod != null ) {
            if( ClassHelper.isDynamicTyped(targetType) || ClassHelper.isDynamicTyped(sourceType) )
                return;
            var lhsOps = resolveOpsType(targetType);
            var rhsOps = resolveOpsType(sourceType);
            var resultType = resolveOpResultType(targetType, sourceType, lhsOps, rhsOps, opMethod);
            if( resultType == null ) {
                addError(String.format("The `%s` operator is not defined for operands with types %s and %s", node.getOperation().getText(), Types.getName(targetType), Types.getName(sourceType)), node);
                return;
            }
            sourceType = resultType;
        }

        if( Types.isAssignableFrom(targetType, sourceType, true) ) {
            if( target instanceof VariableExpression ve && ve.isDynamicTyped() )
                target.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, sourceType);
            else if( target instanceof TupleExpression te )
                applyTupleAssignment(te, sourceType);
        }
        else {
            addError("Assignment target with type " + Types.getName(targetType) + " cannot be assigned to value with type " + Types.getName(sourceType), node);
        }
    }

    /**
     * Get the operator method name for a compound assignment operator
     * (e.g. `+=` corresponds to the `plus` operator), or null if the
     * given operator is not a compound assignment.
     *
     * @param type
     */
    private static String compoundAssignmentMethod(int type) {
        switch( type ) {
            case PLUS_EQUAL:                  return "plus";
            case MINUS_EQUAL:                 return "minus";
            case MULTIPLY_EQUAL:              return "multiply";
            case DIVIDE_EQUAL:                return "div";
            case INTDIV_EQUAL:                return "intdiv";
            case MOD_EQUAL:                   return "mod";
            case POWER_EQUAL:                 return "power";
            case LEFT_SHIFT_EQUAL:            return "leftShift";
            case RIGHT_SHIFT_EQUAL:           return "rightShift";
            case RIGHT_SHIFT_UNSIGNED_EQUAL:  return "rightShiftUnsigned";
            case BITWISE_AND_EQUAL:           return "and";
            case BITWISE_OR_EQUAL:            return "or";
            case BITWISE_XOR_EQUAL:           return "xor";
            default:                          return null;
        }
    }

    private void applyTupleAssignment(TupleExpression target, ClassNode sourceType) {
        if( !TUPLE_TYPE.equals(sourceType) ) {
            addError("Tuple assignment is not compatible with type " + Types.getName(sourceType), target);
            return;
        }

        var vars = target.getExpressions();
        var gts = sourceType.getGenericsTypes();
        if( gts == null )
            return;
        if( vars.size() == gts.length ) {
            for( int i = 0; i < vars.size(); i++ )
                vars.get(i).putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, gts[i].getType());
        }
        else {
            addError("Assignment target with " + vars.size() + " components cannot be assigned to tuple with " + gts.length + " components", target);
        }
    }

    @Override
    public void visitMethodCallExpression(MethodCallExpression node) {
        // resolve argument types (except for closures)
        var arguments = asMethodCallArguments(node);
        boolean hasClosureArgs = false;
        for( var argument : arguments ) {
            if( argument instanceof ClosureExpression ) {
                hasClosureArgs = true;
                continue;
            }
            argument.visit(this);
        }

        // resolve receiver type
        var receiver = node.getObjectExpression();
        if( !node.isImplicitThis() ) {
            receiver.visit(this);
            if( node.isSpreadSafe() ) {
                checkSpreadMethodCall(node);
                return;
            }
            if( ClassHelper.isDynamicTyped(getType(receiver)) ) {
                return;
            }
        }

        // resolve dataflow inputs and outputs for process calls
        if( checkProcessCall(node) )
            return;

        // resolve record return type for record() calls
        if( checkRecordCall(node) )
            return;

        // resolve tuple return type for tuple() calls
        if( checkTupleCall(node) )
            return;

        // resolve method target
        var receiverType = !node.isImplicitThis() ? getType(receiver) : null;
        var target = resolveMethodCall(node);
        if( target != null ) {
            if( arguments.size() > 0 && arguments.get(0) instanceof MapExpression me )
                checkNamedParams(target.getParameters()[0], me);

            if( hasClosureArgs )
                visitClosureArguments(receiverType, arguments, target);

            var conflicts = new ArrayList<GenericsConflict>();
            var dummyMethod = resolveGenericReturnType(receiverType, target, arguments, conflicts);
            checkGenericsConflicts(conflicts);
            node.putNodeMetaData(ASTNodeMarker.METHOD_TARGET, dummyMethod);
            node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, getReturnType(dummyMethod));

            checkOperatorCall(node);
            checkWorkflowCall(node);
        }
        else if( node.getNodeMetaData(ASTNodeMarker.METHOD_TARGET) instanceof MethodNode mn ) {
            var parameters = mn.getParameters();
            if( parameters.length != arguments.size() ) {
                addError(String.format("%s `%s` expects %d argument(s) but received %d", methodType(mn), node.getMethodAsString(), parameters.length, arguments.size()), node.getMethod());
            }
            else {
                checkArguments(parameters, arguments);
                if( arguments.size() > 0 && arguments.get(0) instanceof MapExpression me )
                    checkNamedParams(mn.getParameters()[0], me);
            }

            if( !(mn instanceof ProcessNode || mn instanceof WorkflowNode || mn instanceof AgentNode) ) {
                var dummyMethod = asDummyMethod(receiverType, mn, mn.getParameters(), getReturnType(mn));
                node.putNodeMetaData(ASTNodeMarker.METHOD_TARGET, dummyMethod);
            }
        }
        else if( node.getNodeMetaData(ASTNodeMarker.METHOD_OVERLOADS) != null ) {
            addError(String.format("Function `%s` (with multiple signatures) was called with incorrect number of arguments and/or incorrect argument types", node.getMethodAsString()), node.getMethod());
        }
        else if( !node.isImplicitThis() ) {
            var className = className(receiver);
            addError(String.format("Unrecognized method `%s` for %s", node.getMethodAsString(), className), node.getMethod());
        }
    }

    /**
     * Check a spread-dot method call against an Iterable receiver
     * by resolving the method name and arguments against the
     * receiver's element type.
     *
     * For example, given the following method call:
     *
     *   files('*.txt')*.size()
     *
     * The return type is inferred as follows:
     *
     * 1. The receiver type is Iterable<Path>
     * 2. The element type is Path
     * 3. The target method is Path::size(), which returns an Integer
     * 4. Therefore the return type is Iterable<Integer>
     *
     * @param node
     */
    private void checkSpreadMethodCall(MethodCallExpression node) {
        var receiverType = getType(node.getObjectExpression());
        if( !Types.isAssignableFrom(ClassHelper.ITERABLE_TYPE, receiverType) ) {
            addError("Spread-dot is only supported for Iterable types", node);
            return;
        }
        var elementType = elementType(receiverType);
        var name = node.getMethodAsString();
        var proxyReceiver = new VariableExpression("_", elementType);
        var proxyCall = new MethodCallExpression(proxyReceiver, name, node.getArguments());
        proxyCall.setImplicitThis(false);
        var target = resolveMethodCall(proxyCall);
        if( target != null ) {
            var arguments = asMethodCallArguments(node);
            visitClosureArguments(elementType, arguments, target);

            var conflicts = new ArrayList<GenericsConflict>();
            var dummyMethod = resolveGenericReturnType(elementType, target, arguments, conflicts);
            checkGenericsConflicts(conflicts);
            node.putNodeMetaData(ASTNodeMarker.METHOD_TARGET, dummyMethod);

            var resultType = makeType(receiverType, dummyMethod.getReturnType());
            node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, resultType);
        }
        else {
            addError(String.format("Unrecognized method `%s` for element type %s", name, Types.getName(elementType)), node);
        }
    }

    /**
     * Resolve the return type of dataflow operators where applicable.
     *
     * @param node
     */
    private void checkOperatorCall(MethodCallExpression node) {
        if( node.isImplicitThis() )
            return;

        var receiverType = getType(node.getObjectExpression());
        if( !CHANNEL_TYPE.equals(receiverType) && !VALUE_TYPE.equals(receiverType) )
            return;

        var method = (MethodNode) node.getNodeMetaData(ASTNodeMarker.METHOD_TARGET);
        if( findAnnotation(method, Deprecated.class).isPresent() )
            addSoftError("Operator `" + method.getName() + "` is discouraged from use with static typing", node.getMethod());

        var lhsType = elementType(receiverType);
        var arguments = asMethodCallArguments(node);
        var resultType = new DataflowOpResolver(receiverType, sourceUnit, errorCollector).apply(lhsType, method, arguments);
        if( ClassHelper.isDynamicTyped(resultType) )
            return;

        node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, resultType);
    }

    /**
     * Resolve the return type of a workflow based on the declared outputs.
     *
     * When the workflow declares a single output, the return type contains only
     * that type. When the workflow declares multiple outputs, the return type is
     * a Record containing each named output.
     *
     * @param node
     */
    private void checkWorkflowCall(MethodCallExpression node) {
        var mn = (MethodNode) node.getNodeMetaData(ASTNodeMarker.METHOD_TARGET);
        if( mn instanceof WorkflowNode wn ) {
            var resultType = workflowOutputType(wn.emits);
            node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, resultType);
        }
    }

    private static ClassNode workflowOutputType(Statement block) {
        var emits = asBlockStatements(block);
        if( emits.isEmpty() )
            return ClassHelper.VOID_TYPE;
        if( emits.size() == 1 ) {
            var first = emits.get(0);
            var emit = ((ExpressionStatement) first).getExpression();
            return workflowEmitType(getType(emit));
        }
        var cn = new ClassNode(Record.class);
        for( var stmt : emits ) {
            var emit = ((ExpressionStatement) stmt).getExpression();
            var emitName = outputTarget(emit).getName();
            var emitType = workflowEmitType(getType(emit));
            var fn = new FieldNode(emitName, Modifier.PUBLIC, emitType, cn, null);
            fn.setDeclaringClass(cn);
            cn.addField(fn);
        }
        return cn;
    }

    private static ClassNode workflowEmitType(ClassNode innerType) {
        if( ClassHelper.isDynamicTyped(innerType) )
            return innerType;
        if( CHANNEL_TYPE.equals(innerType) || VALUE_TYPE.equals(innerType) )
            return innerType;
        return makeType(VALUE_TYPE, innerType);
    }

    /**
     * Check the arguments of an invalid method call and report appropriate
     * errors for each invalid argument.
     *
     * @param parameters
     * @param arguments
     */
    private void checkArguments(Parameter[] parameters, List<Expression> arguments) {
        for( int i = 0; i < parameters.length; i++ ) {
            var paramType = parameters[i].getType();
            var argument = arguments.get(i);
            var argType = getType(argument);
            if( Types.isAssignableFrom(paramType, argType) )
                continue;
            if( argument instanceof ClosureExpression ce && Types.isFunctionalInterface(paramType) ) {
                var parameterTypes = Arrays.stream(ce.getParameters())
                    .map(p -> p.getType())
                    .toArray(ClassNode[]::new);
                var returnType = ClassHelper.dynamicType();
                addError("Closure with signature " + Types.getName(parameterTypes, returnType) + " is not compatible with expected signature: " + Types.getName(paramType), argument);
            }
            else {
                addArgumentTypeError(argument, argType, paramType);
            }
        }
    }

    /**
     * Report call arguments that inferred conflicting types for the same
     * type parameter. For example, `channel.of(1, 'a')` infers the element
     * type as both Integer and String.
     *
     * @param conflicts
     */
    private void checkGenericsConflicts(List<GenericsConflict> conflicts) {
        for( var conflict : conflicts )
            addArgumentTypeError(conflict.argument(), conflict.actualType(), conflict.expectedType());
    }

    private void addArgumentTypeError(Expression argument, ClassNode argType, ClassNode paramType) {
        addError("Argument with type " + Types.getName(argType) + " is not compatible with parameter of type " + Types.getName(paramType), argument);
    }

    /**
     * Check named arguments against the named parameters of a
     * method, if they are defined using the @NamedParams annotation.
     *
     * @param param
     * @param args
     */
    private void checkNamedParams(Parameter param, MapExpression args) {
        // compare raw types -- the named-params map carries injected field info as
        // type arguments, and MAP_TYPE carries its K,V placeholders
        if( !Types.isEqual(param.getType().getPlainNodeReference(), ClassHelper.MAP_TYPE.getPlainNodeReference()) )
            return;

        var namedParams = asNamedParams(param);
        if( namedParams.isEmpty() )
            return;

        for( var entry : args.getMapEntryExpressions() ) {
            var name = entry.getKeyExpression().getText();
            var value = entry.getValueExpression();
            if( !namedParams.containsKey(name) ) {
                addError("Named param `" + name + "` is not defined", entry);
                continue;
            }
            var namedParam = asNamedParam(namedParams.get(name));
            var argType = getType(value);
            if( !Types.isAssignableFrom(namedParam.getType(), argType) )
                addError("Named param `" + name + "` expects a " + Types.getName(namedParam.getType()) + " but received a " + Types.getName(argType), value);
            entry.putNodeMetaData("_NAMED_PARAM", namedParam);
        }
    }

    /**
     * Visit closure arguments after the target method is resolved,
     * so that closure parameters can be inferred from the corresponding
     * method parameter.
     *
     * @param receiverType
     * @param arguments
     * @param method
     */
    private void visitClosureArguments(ClassNode receiverType, List<Expression> arguments, MethodNode method) {
        var parameters = method.getParameters();
        for( int i = 0; i < arguments.size(); i++ ) {
            var argument = arguments.get(i);
            if( !(argument instanceof ClosureExpression source) )
                continue;
            var target = parameters[i];
            var resolvedPlaceholders = resolveGenericsConnections(receiverType, method, arguments);
            var samParameterTypes = resolveClosureParameterTypes(source, target, resolvedPlaceholders);
            argument.visit(this);
        }
    }

    /**
     * Check process calls against the declared process inputs, and resolve
     * the return type based on the declared process outputs.
     *
     * When calling a process, each argument can be either the same type as
     * the declared input, channel of that type, or a dataflow value of that
     * type. For example, given a declared input with type String, the argument
     * should be either a String, Channel<String>, or Value<String>.
     *
     * The return type is determined by the argument types and the declared outputs.
     *
     * - When a process is called with a Channel argument, the output is wrapped
     *   in a Channel, otherwise it is wrapped in a Value.
     *
     * - When a process declares a single output, the return type contains only that
     *   type (e.g. T -> Channel<T>). When a process declares multiple outputs,
     *   the return type is a Record where each field is wrapped in the appropriate
     *   dataflow type.
     *
     * @param node
     */
    private boolean checkProcessCall(MethodCallExpression node) {
        var mn = (MethodNode) node.getNodeMetaData(ASTNodeMarker.METHOD_TARGET);
        var outputs = processOutputs(mn);
        if( outputs == null )
            return false;
        var label = methodType(mn);

        var parameters = mn.getParameters();
        var arguments = asMethodCallArguments(node);
        if( parameters.length != arguments.size() )
            return false;

        for( int i = 0; i < parameters.length; i++ ) {
            var paramType = parameters[i].getType();
            var argType = getType(arguments.get(i));
            var elementType = dataflowElementType(argType);
            if( !Types.isAssignableFrom(paramType, elementType) )
                addError("Argument with type " + Types.getName(elementType) + " is not compatible with " + label.toLowerCase() + " input of type " + Types.getName(paramType), arguments.get(i));
        }

        var numChannelArgs = arguments.stream().filter((arg) -> CHANNEL_TYPE.equals(getType(arg))).count();
        var numDynamicArgs = arguments.stream().filter((arg) -> ClassHelper.isDynamicTyped(getType(arg))).count();
        if( numChannelArgs > 1 )
            addError(label + " `" + mn.getName() + "` was called with multiple channel arguments which can lead to non-deterministic behavior -- make sure that at most one argument is a channel and that all other arguments are dataflow values", node);

        var dataflowType = numChannelArgs + numDynamicArgs > 0 ? CHANNEL_TYPE : VALUE_TYPE;
        var resultType = processOutputType(dataflowType, outputs);
        node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, resultType);

        return true;
    }

    /**
     * The output block of a process or agent, or null for anything else.
     * An agent has the same typed inputs and outputs as a typed process,
     * so it gets the same dataflow treatment at the call site.
     */
    private static Statement processOutputs(MethodNode mn) {
        if( mn instanceof ProcessNodeV2 pn )
            return pn.outputs;
        if( mn instanceof AgentNode an )
            return an.outputs;
        return null;
    }

    private static ClassNode dataflowElementType(ClassNode type) {
        if( CHANNEL_TYPE.equals(type) || VALUE_TYPE.equals(type) )
            return elementType(type);
        return type;
    }

    private static ClassNode processOutputType(ClassNode dataflowType, Statement block) {
        var outputs = asBlockStatements(block);
        if( outputs.isEmpty() )
            return ClassHelper.VOID_TYPE;
        if( outputs.size() == 1 ) {
            var first = outputs.get(0);
            var output = ((ExpressionStatement) first).getExpression();
            return processEmitType(dataflowType, getType(output));
        }
        var cn = new ClassNode(Record.class);
        for( var stmt : outputs ) {
            var output = ((ExpressionStatement) stmt).getExpression();
            var emitName = outputTarget(output).getName();
            var emitType = processEmitType(dataflowType, getType(output));
            var fn = new FieldNode(emitName, Modifier.PUBLIC, emitType, cn, null);
            fn.setDeclaringClass(cn);
            cn.addField(fn);
        }
        return cn;
    }

    private static ClassNode processEmitType(ClassNode dataflowType, ClassNode innerType) {
        return makeType(dataflowType, innerType);
    }

    private static VariableExpression outputTarget(Expression output) {
        if( output instanceof VariableExpression ve )
            return ve;
        if( output instanceof AssignmentExpression ae )
            return (VariableExpression)ae.getLeftExpression();
        return null;
    }

    private static String methodType(MethodNode node) {
        if( node instanceof ProcessNode )
            return "Process";
        if( node instanceof AgentNode )
            return "Agent";
        if( node instanceof WorkflowNode )
            return "Workflow";
        return "Function";
    }

    private static String className(Expression node) {
        var receiverType = getType(node);
        return receiverType != null && receiverType.implementsInterface(ClassHelper.makeCached(Namespace.class))
            ? "namespace `" + node.getText() + "`"
            : "type " + Types.getName(receiverType);
    }

    /**
     * Resolve record() calls by creating a custom Record type
     * with fields corresponding to the call arguments.
     *
     * For example, the return type of `record(id: '1', count: 42)`
     * should be Record<id: String, count: Integer>.
     *
     * @param node
     */
    private boolean checkRecordCall(MethodCallExpression node) {
        if( !node.isImplicitThis() || !"record".equals(node.getMethodAsString()) )
            return false;

        var mn = (MethodNode) node.getNodeMetaData(ASTNodeMarker.METHOD_TARGET);
        if( mn == null || !RECORD_TYPE.equals(mn.getReturnType()) )
            return false;

        var returnType = new ClassNode(Record.class);
        for( var entry : asNamedArgs(node) ) {
            var name = entry.getKeyExpression().getText();
            var type = getType(entry.getValueExpression());
            var fn = new FieldNode(name, Modifier.PUBLIC, type, returnType, null);
            fn.setDeclaringClass(returnType);
            returnType.addField(fn);
        }

        node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, returnType);
        return true;
    }

    /**
     * Resolve tuple() calls by creating a custom Tuple type
     * that is parameterized based on the call arguments.
     *
     * For example, the return type of `tuple('hello', 42, true)`
     * should be Tuple<String, Integer, Boolean>.
     *
     * @param node
     */
    private boolean checkTupleCall(MethodCallExpression node) {
        if( !node.isImplicitThis() || !"tuple".equals(node.getMethodAsString()) )
            return false;

        var mn = (MethodNode) node.getNodeMetaData(ASTNodeMarker.METHOD_TARGET);
        if( mn == null || !TUPLE_TYPE.equals(mn.getReturnType()) )
            return false;

        var arguments = asMethodCallArguments(node);
        var dummyParams = IntStream.range(0, arguments.size())
            .mapToObj(i -> new Parameter(getType(arguments.get(i)), "arg" + i))
            .toArray(Parameter[]::new);
        var returnType = new ClassNode(Tuple.class);
        var genericsTypes = arguments.stream()
            .map(arg -> new GenericsType(getType(arg)))
            .toArray(GenericsType[]::new);
        returnType.setGenericsTypes(genericsTypes);

        var dummyMethod = asDummyMethod(null, mn, dummyParams, returnType);
        node.putNodeMetaData(ASTNodeMarker.METHOD_TARGET, dummyMethod);
        node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, returnType);
        return true;
    }

    @Override
    public void visitConstructorCallExpression(ConstructorCallExpression node) {
        var type = node.getType();
        if( Types.isRecordType(type) )
            addError("Record type " + Types.getName(type) + " cannot be used as a constructor -- use `record()` instead", node);
        super.visitConstructorCallExpression(node);
    }

    @Override
    public void visitBinaryExpression(BinaryExpression node) {
        super.visitBinaryExpression(node);

        var op = node.getOperation();
        if( op.getType() == PLUS && checkRecordSum(node) )
            return;
        if( op.getType() == LEFT_SQUARE_BRACKET && checkTupleComponent(node) )
            return;

        var lhsType = getType(node.getLeftExpression());
        var rhsType = getType(node.getRightExpression());
        if( ClassHelper.isDynamicTyped(lhsType) || ClassHelper.isDynamicTyped(rhsType) )
            return;

        var lhsOps = resolveOpsType(lhsType);
        var rhsOps = resolveOpsType(rhsType);

        ClassNode resultType = null;
        switch( op.getType() ) {
            case POWER:
                resultType = resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "power");
                break;

            case MULTIPLY:
                resultType = resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "multiply");
                break;

            case DIVIDE:
                resultType = resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "div");
                break;

            case MOD:
                resultType = resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "mod");
                break;

            case PLUS:
                resultType = resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "plus");
                break;

            case MINUS:
                resultType = resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "minus");
                break;

            case LEFT_SHIFT:
                resultType = resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "leftShift");
                break;

            case RIGHT_SHIFT:
                resultType = resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "rightShift");
                break;

            case RIGHT_SHIFT_UNSIGNED:
                resultType = resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "rightShiftUnsigned");
                break;

            case KEYWORD_INSTANCEOF:
            case COMPARE_NOT_INSTANCEOF:
                return;

            case COMPARE_TO:
                resultType = resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "compareTo");
                break;

            case COMPARE_LESS_THAN:
            case COMPARE_LESS_THAN_EQUAL:
            case COMPARE_GREATER_THAN:
            case COMPARE_GREATER_THAN_EQUAL:
                resultType = resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "compareTo") != null ? ClassHelper.Boolean_TYPE : null;
                break;

            case KEYWORD_IN:
            case COMPARE_NOT_IN:
                resultType = resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "isCase");
                break;

            case COMPARE_EQUAL:
            case COMPARE_NOT_EQUAL:
                resultType = Types.isEqual(lhsType, rhsType) || resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "compareTo") != null
                    ? ClassHelper.Boolean_TYPE
                    : null;
                break;

            case FIND_REGEX:
                resultType = Types.isEqual(ClassHelper.STRING_TYPE, lhsType) && Types.isEqual(ClassHelper.STRING_TYPE, rhsType) ? ClassHelper.dynamicType() : null;
                break;

            case MATCH_REGEX:
                resultType = Types.isEqual(ClassHelper.STRING_TYPE, lhsType) && Types.isEqual(ClassHelper.STRING_TYPE, rhsType) ? ClassHelper.Boolean_TYPE : null;
                break;

            case BITWISE_AND:
                resultType = resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "and");
                break;

            case BITWISE_XOR:
                resultType = resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "xor");
                break;

            case BITWISE_OR:
                resultType = resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "or");
                break;

            case LOGICAL_AND:
            case LOGICAL_OR:
                return;

            case LEFT_SQUARE_BRACKET:
                resultType = resolveOpResultType(lhsType, rhsType, lhsOps, rhsOps, "getAt");
                break;
        }

        if( resultType != null )
            node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, resultType);
        else
            addError(String.format("The `%s` operator is not defined for operands with types %s and %s", op.getText(), Types.getName(lhsType), Types.getName(rhsType)), node);
    }

    /**
     * Resolve the type of a record sum.
     *
     * For example, the expression `record(foo: 1) + record(bar: '...')` has type
     * Record { foo: Integer ; bar: String }.
     *
     * @param node
     */
    private boolean checkRecordSum(BinaryExpression node) {
        var lhs = node.getLeftExpression();
        var rhs = node.getRightExpression();
        var lhsType = getType(lhs);
        var rhsType = getType(rhs);
        if( !Types.isRecordType(lhsType) || !Types.isRecordType(rhsType) )
            return false;

        var resultType = recordSumType(lhsType, rhsType);
        node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, resultType);
        return true;
    }

    /**
     * Resolve the type of a tuple component expression.
     *
     * For example, the expression `tuple('hello', 42, true)[1]` has type Integer.
     *
     * @param node
     */
    private boolean checkTupleComponent(BinaryExpression node) {
        var lhs = node.getLeftExpression();
        var rhs = node.getRightExpression();
        var lhsType = getType(lhs);
        if( !TUPLE_TYPE.equals(lhsType) )
            return false;
        var index = rhs instanceof ConstantExpression ce && ce.getValue() instanceof Integer i ? i : null;
        if( index == null ) {
            addError("Tuple component index should be an integer literal", rhs);
            return true;
        }
        var gts = lhsType.getGenericsTypes();
        if( gts == null )
            return true;
        if( index < 0 || index >= gts.length ) {
            addError("Tuple component index is out of range -- it should range between 0 and " + (gts.length - 1), node);
            return true;
        }
        node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, gts[index].getType());
        return true;
    }

    @Override
    public void visitTernaryExpression(TernaryExpression node) {
        super.visitTernaryExpression(node);
        applyConditionalExpression(node);
    }

    @Override
    public void visitShortTernaryExpression(ElvisOperatorExpression node) {
        super.visitShortTernaryExpression(node);
        applyConditionalExpression(node);
    }

    private void applyConditionalExpression(TernaryExpression node) {
        var trueExpr = node.getTrueExpression();
        var falseExpr = node.getFalseExpression();
        var trueType = getType(trueExpr);
        var falseType = getType(falseExpr);

        ClassNode resultType;
        boolean nullable = true;
        if( isNullConstant(trueExpr) && isNullConstant(falseExpr) ) {
            resultType = null;
        }
        else if( !isNullConstant(trueExpr) && isNullConstant(falseExpr) ) {
            resultType = trueType;
        }
        else if( isNullConstant(trueExpr) && !isNullConstant(falseExpr) ) {
            resultType = falseType;
        }
        else if( ClassHelper.isDynamicTyped(trueType) || ClassHelper.isDynamicTyped(falseType) ) {
            return;
        }
        else if( Types.isEqual(trueType, falseType) ) {
            resultType = trueType;
            nullable = isNullable(trueType) || isNullable(falseType);
        }
        else if( isEmptyListOrMap(falseExpr) && Types.isAssignableFrom(trueType, falseType) ) {
            resultType = trueType;
            nullable = isNullable(trueType);
        }
        else if( isEmptyListOrMap(trueExpr) && Types.isAssignableFrom(falseType, trueType) ) {
            resultType = falseType;
            nullable = isNullable(falseType);
        }
        else {
            addError(String.format("Conditional expression has inconsistent types -- true branch has type %s but false branch has type %s", Types.getName(trueType), Types.getName(falseType)), node);
            return;
        }

        if( nullable && resultType != null ) {
            resultType = new ClassNode(resultType.getTypeClass());
            resultType.putNodeMetaData(ASTNodeMarker.NULLABLE, Boolean.TRUE);
        }
        node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, resultType);
    }

    private static boolean isEmptyListOrMap(Expression node) {
        if( node instanceof ListExpression le )
            return le.getExpressions().isEmpty();
        if( node instanceof MapExpression me )
            return me.getMapEntryExpressions().isEmpty();
        return false;
    }

    private static boolean isNullable(ClassNode cn) {
        return cn == null || cn.getNodeMetaData(ASTNodeMarker.NULLABLE) != null;
    }

    @Override
    public void visitClosureExpression(ClosureExpression node) {
        super.visitClosureExpression(node);

        // resolve return type and check against declared return type. A closure
        // always returns its last expression, and a void signature discards it,
        // so there is nothing to check against. Likewise, any value can be
        // coerced to a boolean, so a predicate can return any type.
        var returnType = (ClassNode) node.getNodeMetaData(ASTNodeMarker.INFERRED_RETURN_TYPE);
        if( returnType != null && !returnType.equals(ClassHelper.VOID_TYPE) && !Types.isEqual(returnType, ClassHelper.Boolean_TYPE) ) {
            var visitor = new ReturnStatementVisitor(sourceUnit, errorCollector);
            visitor.visit(returnType, node.getCode());

            var inferredReturnType = visitor.getInferredReturnType();
            if( inferredReturnType != null )
                node.putNodeMetaData(ASTNodeMarker.INFERRED_RETURN_TYPE, inferredReturnType);
        }
    }

    @Override
    public void visitListExpression(ListExpression node) {
        super.visitListExpression(node);

        ClassNode elementType = null;
        for( var el : node.getExpressions() ) {
            var type = getType(el);
            if( elementType == null ) {
                elementType = type;
            }
            else if( !Types.isEqual(elementType, type) ) {
                addError(String.format("List expression has inconsistent element types -- some elements have type %s while others have type %s", Types.getName(elementType), Types.getName(type)), node);
                break;
            }
        }

        if( elementType != null ) {
            var resultType = makeType(ClassHelper.LIST_TYPE, elementType);
            node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, resultType);
        }
    }

    @Override
    public void visitMapExpression(MapExpression node) {
        super.visitMapExpression(node);

        var entries = node.getMapEntryExpressions();
        if( entries.isEmpty() ) {
            node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, ClassHelper.MAP_TYPE.getPlainNodeReference());
            return;
        }

        // infer key and value type args, falling back to `?` when heterogeneous
        var keyType = commonType(entries.stream().map(e -> getType(e.getKeyExpression())));
        var valueType = commonType(entries.stream().map(e -> getType(e.getValueExpression())));
        var resultType = makeType(ClassHelper.MAP_TYPE, keyType, valueType);
        node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, resultType);
    }

    private static ClassNode commonType(Stream<ClassNode> types) {
        return types.reduce((a, b) -> Types.isEqual(a, b) ? a : ClassHelper.dynamicType())
            .orElse(ClassHelper.dynamicType());
    }

    @Override
    public void visitRangeExpression(RangeExpression node) {
        super.visitRangeExpression(node);

        var lhs = node.getFrom();
        var rhs = node.getTo();
        var lhsType = getType(lhs);
        var rhsType = getType(rhs);

        if( !Types.isEqual(lhsType, rhsType) ) {
            addError("Lower bound and upper bound of range expression should have the same type", node);
            node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, ClassHelper.dynamicType());
            return;
        }

        var elementType = Arrays.stream(RANGE_TYPES)
            .filter(type -> Types.isEqual(type, lhsType))
            .findFirst().orElse(null);

        if( elementType != null ) {
            var resultType = makeType(ClassHelper.LIST_TYPE, elementType);
            node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, resultType);
        }
        else {
            addError("Range expression with elements of type " + Types.getName(lhsType) + " is not supported", node);
            node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, ClassHelper.dynamicType());
        }
    }

    private static final ClassNode[] RANGE_TYPES = new ClassNode[] { ClassHelper.Integer_TYPE, ClassHelper.STRING_TYPE };

    @Override
    public void visitUnaryMinusExpression(UnaryMinusExpression node) {
        super.visitUnaryMinusExpression(node);
        resolveUnaryOpOrFail(node.getExpression(), "-", "negative", node);
    }

    @Override
    public void visitUnaryPlusExpression(UnaryPlusExpression node) {
        super.visitUnaryPlusExpression(node);
        resolveUnaryOpOrFail(node.getExpression(), "+", "positive", node);
    }

    @Override
    public void visitBitwiseNegationExpression(BitwiseNegationExpression node) {
        super.visitBitwiseNegationExpression(node);
        resolveUnaryOpOrFail(node.getExpression(), "~", "bitwiseNegate", node);
    }

    private void resolveUnaryOpOrFail(Expression operand, String op, String method, ASTNode node) {
        var type = getType(operand);
        var opsType = resolveOpsType(type);
        var resultType = resolveOpResultType(opsType, method);
        if( resultType != null )
            node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, resultType);
        else
            addError(String.format("The `%s` operator is not defined for an operand with type %s", op, Types.getName(type)), node);
    }

    @Override
    public void visitCastExpression(CastExpression node) {
        super.visitCastExpression(node);
        var sourceType = getType(node.getExpression());
        var targetType = node.getType();
        if( ClassHelper.isObjectType(sourceType) || Types.isAssignableFrom(targetType, sourceType) )
            return;
        if( checkRecordCast(targetType, sourceType, node) )
            return;
        var opsType = resolveOpsType(targetType);
        if( resolveOpResultType(sourceType, opsType, "ofType") == null )
            addError(String.format("Value of type %s cannot be cast to %s", Types.getName(sourceType), Types.getName(targetType)), node);
    }

    private boolean checkRecordCast(ClassNode targetType, ClassNode sourceType, ASTNode node) {
        if( !(targetType.redirect() instanceof RecordNode) )
            return false;
        if( !Types.isRecordType(sourceType) )
            return false;
        for( var target : targetType.getFields() ) {
            if( target.getType().getNodeMetaData(ASTNodeMarker.NULLABLE) != null )
                continue;
            var source = sourceType.getDeclaredField(target.getName());
            if( source == null ) {
                addError(String.format("Record mismatch -- source record is missing field `%s` required by %s", target.getName(), Types.getName(targetType)), node);
                continue;
            }
            if( !Types.isAssignableFrom(target.getType(), source.getType()) ) {
                addError(String.format("Record mismatch -- field `%s` of %s expects a %s but received a %s", target.getName(), Types.getName(targetType), Types.getName(target.getType()), Types.getName(source.getType())), node);
                continue;
            }
        }
        return true;
    }

    @Override
    public void visitPropertyExpression(PropertyExpression node) {
        super.visitPropertyExpression(node);

        var mn = asMethodOutput(node);
        if( mn instanceof ProcessNode || mn instanceof WorkflowNode ) {
            addError("Using the `.out` property to access process/workflow outputs is not supported with static typing -- assign the output to a variable instead", node);
            return;
        }

        var receiver = node.getObjectExpression();
        var receiverType = getType(receiver);
        if( ClassHelper.isDynamicTyped(receiverType) )
            return;
        if( RECORD_TYPE.equals(receiverType) && receiverType.getFields().isEmpty() )
            return;

        if( node.isSpreadSafe() ) {
            checkSpreadProperty(node);
            return;
        }

        var target = resolveProperty(node);
        if( target != null ) {
            node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, target.getType());
            return;
        }

        var property = node.getPropertyAsString();
        if( Types.isEqual(receiverType, PARAMS_TYPE) ) {
            if( hasParamsBlock )
                addError("Unrecognized parameter `" + property + "`", node);
        }
        else {
            var className = className(receiver);
            addError(String.format("Unrecognized property `%s` for %s", property, className), node);
        }
    }

    /**
     * Check a spread-dot property access against an Iterable receiver
     * by resolving the method name and arguments against the
     * receiver's element type.
     *
     * For example, given the following property access:
     *
     *   files('*.txt')*.name
     *
     * The result type is inferred as follows:
     *
     * 1. The receiver type is Iterable<Path>
     * 2. The element type is Path
     * 3. The target field is Path::name, whose type is String
     * 4. Therefore the result type is Iterable<String>
     *
     * @param node
     */
    private void checkSpreadProperty(PropertyExpression node) {
        var receiverType = getType(node.getObjectExpression());
        if( !Types.isAssignableFrom(ClassHelper.ITERABLE_TYPE, receiverType) ) {
            addError("Spread-dot is only supported for Iterable types", node);
            return;
        }
        var elementType = elementType(receiverType);
        if( RECORD_TYPE.equals(elementType) && elementType.getFields().isEmpty() )
            return;
        var property = node.getPropertyAsString();
        var fn = resolveProperty(elementType, property);
        if( fn != null ) {
            var resultType = makeType(receiverType, fn.getType());
            node.putNodeMetaData(ASTNodeMarker.PROPERTY_TARGET, fn);
            node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, resultType);
        }
        else {
            addError(String.format("Unrecognized property `%s` for element type %s", property, Types.getName(elementType)), node);
        }
    }

    public void addWarning(String message, String tokenText, ASTNode node) {
        var token = new Token(0, tokenText, node.getLineNumber(), node.getColumnNumber()); // ASTNode to CSTNode
        errorCollector.addWarning(WarningMessage.POSSIBLE_ERRORS, message, token, sourceUnit);
    }

    @Override
    public void addError(String message, ASTNode node) {
        addError(new TypeError(message, node));
    }

    /**
     * Report an error that should not fail a lint. A legacy operator still works,
     * so flagging one is advice, not a defect.
     */
    public void addSoftError(String message, ASTNode node) {
        addError(new TypeError(message, node, true));
    }

    private void addError(TypeError cause) {
        var errorMessage = new SyntaxErrorMessage(cause, sourceUnit);
        errorCollector.addErrorAndContinue(errorMessage);
    }

}
