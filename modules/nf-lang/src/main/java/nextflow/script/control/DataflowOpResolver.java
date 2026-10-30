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
package nextflow.script.control;

import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.types.Bag;
import nextflow.script.types.Channel;
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
import org.codehaus.groovy.ast.expr.BinaryExpression;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.ConstantExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.NamedArgumentListExpression;
import org.codehaus.groovy.ast.expr.PropertyExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.ReturnStatement;
import org.codehaus.groovy.control.ErrorCollector;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;

import static nextflow.script.ast.ASTUtils.*;
import static nextflow.script.types.TypeCheckingUtils.*;

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class DataflowOpResolver {

    private static final ClassNode BAG_TYPE = ClassHelper.makeCached(Bag.class);
    private static final ClassNode CHANNEL_TYPE = ClassHelper.makeCached(Channel.class);
    private static final ClassNode TUPLE_TYPE = ClassHelper.makeCached(Tuple.class);
    private static final ClassNode VALUE_TYPE = ClassHelper.makeCached(Value.class);

    private ClassNode receiverType;
    private SourceUnit sourceUnit;
    private ErrorCollector errorCollector;

    public DataflowOpResolver(ClassNode receiverType, SourceUnit sourceUnit, ErrorCollector errorCollector) {
        this.receiverType = receiverType;
        this.sourceUnit = sourceUnit;
        this.errorCollector = errorCollector;
    }

    /**
     * Resolve the return type of dataflow operators where applicable,
     * such as `combine`, `groupBy`, and `join`.
     *
     * @param lhsType
     * @param method
     * @param arguments
     */
    public ClassNode apply(ClassNode lhsType, MethodNode method, List<Expression> arguments) {
        var name = method.getName();

        if( "combine".equals(name) )
            return applyCombine(lhsType, arguments);

        if( "groupBy".equals(name) )
            return applyGroupBy(lhsType, arguments);

        if( "join".equals(name) )
            return applyJoin(lhsType, arguments);

        return ClassHelper.dynamicType();
    }

    /**
     * Infer the parameter types of a closure given as a named argument
     * to a dataflow operator, before the closure is visited.
     *
     * @param lhsType
     * @param method
     * @param name
     * @param closure
     * @param arguments
     */
    public void inferNamedClosureParameterTypes(ClassNode lhsType, MethodNode method, String name, ClosureExpression closure, List<Expression> arguments) {
        if( "join".equals(method.getName()) && "by".equals(name) )
            inferJoinClosureParameterTypes(lhsType, closure, arguments);
    }

    /**
     * Resolve the result type of a `combine` operation in terms of the left
     * and right operands.
     *
     * Given arguments of type `L` and `R`, `combine` produces a tuple of type
     * `(L, R)`. If `L` and/or `R` are tuples, they are flattened into the resulting
     * tuple.
     *
     * @param lhsType
     * @param arguments
     */
    private ClassNode applyCombine(ClassNode lhsType, List<Expression> arguments) {
        if( arguments.size() == 1 && arguments.get(0) instanceof NamedArgumentListExpression nale )
            return applyCombineNamedArgs(lhsType, nale);
        var rhsType = dataflowElementType(getType(arguments.get(0)));
        var componentTypes = new ArrayList<ClassNode>();
        if( !combineTupleOrValue(componentTypes, lhsType) )
            return channelTupleType(null);
        if( !combineTupleOrValue(componentTypes, rhsType) )
            return channelTupleType(null);
        var gts = componentTypes.stream()
            .map(cn -> new GenericsType(cn))
            .toArray(GenericsType[]::new);
        return channelTupleType(gts);
    }

    private ClassNode applyCombineNamedArgs(ClassNode lhsType, NamedArgumentListExpression nale) {
        if( !Types.isRecordType(lhsType) )
            return ClassHelper.dynamicType();
        var rhsType = new ClassNode(Record.class);
        for( var entry : nale.getMapEntryExpressions() ) {
            var name = entry.getKeyExpression().getText();
            var value = entry.getValueExpression();
            var valueType = dataflowValueType(getType(value), value);
            var fn = new FieldNode(name, Modifier.PUBLIC, valueType, rhsType, null);
            fn.setDeclaringClass(rhsType);
            rhsType.addField(fn);
        }
        var elementType = recordSumType(lhsType, rhsType);
        return makeType(receiverType, elementType);
    }

    private ClassNode dataflowValueType(ClassNode type, ASTNode node) {
        if( CHANNEL_TYPE.equals(type) ) {
            addError("Named argument of a `combine` operation cannot be a channel -- it must be a dataflow value or plain value", node);
            return ClassHelper.dynamicType();
        }
        if( VALUE_TYPE.equals(type) ) {
            return elementType(type);
        }
        return type;
    }

    private static boolean combineTupleOrValue(List<ClassNode> componentTypes, ClassNode type) {
        if( TUPLE_TYPE.equals(type) ) {
            var gts = type.getGenericsTypes();
            if( gts == null || gts.length == 0 )
                return false;
            for( int i = 0; i < gts.length; i++ )
                componentTypes.add(gts[i].getType());
        }
        else {
            componentTypes.add(type);
        }
        return true;
    }

    private ClassNode channelTupleType(GenericsType[] gts) {
        var tupleType = TUPLE_TYPE.getPlainNodeReference();
        tupleType.setGenericsTypes(gts);
        return makeType(receiverType, tupleType);
    }

    /**
     * Resolve the result type of a `groupBy` operation.
     *
     * Given source tuples of type `(K, N, V)` or `(K, V)`,
     * `groupBy` produces a tuple of type `(K, Bag<V>)`.
     *
     * @param lhsType
     * @param arguments
     */
    private ClassNode applyGroupBy(ClassNode lhsType, List<Expression> arguments) {
        if( !TUPLE_TYPE.equals(lhsType) )
            return ClassHelper.dynamicType();
        var lgts = lhsType.getGenericsTypes();
        if( lgts == null || !(lgts.length == 2 || lgts.length == 3) )
            return ClassHelper.dynamicType();
        var keyType = lgts[0].getType();
        var valueType = lgts[lgts.length - 1].getType();
        var gts = new GenericsType[] {
            new GenericsType(keyType),
            new GenericsType(makeType(BAG_TYPE, valueType))
        };
        return channelTupleType(gts);
    }

    /**
     * Resolve the result type of a `join` operation in terms of the left
     * and right operands.
     *
     * Given two matching records R1 and R2, `join` produces R1 + R2.
     *
     * @param lhsType
     * @param arguments
     */
    private ClassNode applyJoin(ClassNode lhsType, List<Expression> arguments) {
        if( !Types.isRecordType(lhsType) )
            return ClassHelper.dynamicType();
        var argType = getType(arguments.get(arguments.size() - 1));
        var rhsType = dataflowElementType(argType);
        if( !Types.isRecordType(rhsType) )
            return ClassHelper.dynamicType();
        checkJoinFields(lhsType, rhsType, arguments);
        var elementType = recordSumType(lhsType, rhsType);
        return makeType(CHANNEL_TYPE, elementType);
    }

    /**
     * Infer the parameter types of a `by` closure for a `join` operation.
     *
     * A key function is applied to records from both sides, so its parameter
     * is typed as the fields common to both records. A join predicate receives
     * a record from each side.
     *
     * @param lhsType
     * @param closure
     * @param arguments
     */
    private void inferJoinClosureParameterTypes(ClassNode lhsType, ClosureExpression closure, List<Expression> arguments) {
        var rhsType = dataflowElementType(getType(arguments.get(arguments.size() - 1)));
        if( !Types.isRecordType(lhsType) || !Types.isRecordType(rhsType) )
            return;
        var params = closure.getParameters();
        if( params == null )
            return;
        if( params.length == 1 ) {
            applyParameterType(params[0], recordCommonType(lhsType, rhsType));
        }
        else if( params.length == 2 ) {
            applyParameterType(params[0], lhsType);
            applyParameterType(params[1], rhsType);
            closure.putNodeMetaData(ASTNodeMarker.INFERRED_RETURN_TYPE, ClassHelper.Boolean_TYPE);
        }
    }

    private static void applyParameterType(Parameter parameter, ClassNode type) {
        if( parameter.isDynamicTyped() )
            parameter.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, type);
    }

    /**
     * Get the record type containing the fields that are present
     * in both records. Nested records are intersected recursively,
     * and fields with different types are dynamically typed.
     *
     * @param lhsType
     * @param rhsType
     */
    private static ClassNode recordCommonType(ClassNode lhsType, ClassNode rhsType) {
        var resultType = new ClassNode(Record.class);
        for( var lfn : lhsType.getFields() ) {
            var rfn = rhsType.getField(lfn.getName());
            if( rfn == null )
                continue;
            var fieldType = Types.isRecordType(lfn.getType()) && Types.isRecordType(rfn.getType())
                ? recordCommonType(lfn.getType(), rfn.getType())
                : Types.isEqual(lfn.getType(), rfn.getType()) ? lfn.getType() : ClassHelper.dynamicType();
            var fn = new FieldNode(lfn.getName(), Modifier.PUBLIC, fieldType, resultType, null);
            fn.setDeclaringClass(resultType);
            resultType.addField(fn);
        }
        return resultType;
    }

    private void checkJoinFields(ClassNode lhsType, ClassNode rhsType, List<Expression> arguments) {
        if( arguments.isEmpty() || !(arguments.get(0) instanceof NamedArgumentListExpression nale) )
            return;
        for( var entry : nale.getMapEntryExpressions() ) {
            if( !"by".equals(entry.getKeyExpression().getText()) )
                continue;
            var value = entry.getValueExpression();
            if( value instanceof ClosureExpression closure )
                checkJoinClosure(lhsType, rhsType, closure);
            else if( value instanceof ConstantExpression ce && ce.getValue() instanceof String name )
                checkJoinField(lhsType, rhsType, name, value);
            else if( !ClassHelper.isDynamicTyped(getType(value)) && !Types.isEqual(getType(value), ClassHelper.STRING_TYPE) )
                addError("Join option `by` must be a field name or a closure", value);
        }
    }

    private void checkJoinField(ClassNode lhsType, ClassNode rhsType, String name, Expression field) {
        // Report an error if the `by` field is not present in both records
        var lhsMissing = lhsType.getField(name) == null;
        var rhsMissing = rhsType.getField(name) == null;
        if( lhsMissing )
            addError("Join field `" + name + "` is not present in left-hand side", field);
        if( rhsMissing )
            addError("Join field `" + name + "` is not present in right-hand side", field);
        if( lhsMissing || rhsMissing )
            return;
        checkOverwrittenFields(lhsType, rhsType, Set.of(name), field);
    }

    private void checkJoinClosure(ClassNode lhsType, ClassNode rhsType, ClosureExpression closure) {
        var params = closure.getParameters();
        var arity = params == null ? 0 : closure.isParameterSpecified() ? params.length : 1;
        if( arity == 1 )
            checkOverwrittenFields(lhsType, rhsType, keyFunctionFields(closure), closure);
        else if( arity == 2 )
            checkOverwrittenFields(lhsType, rhsType, joinPredicateFields(closure), closure);
        else
            addError("Join closure must have one parameter (key function) or two parameters (join predicate)", closure);
    }

    private void checkOverwrittenFields(ClassNode lhsType, ClassNode rhsType, Set<String> keyFields, ASTNode node) {
        // Report a soft error if a field other than the matching key is present in both records
        var overwritten = rhsType.getFields().stream()
            .map(fn -> fn.getName())
            .filter(fieldName -> !keyFields.contains(fieldName) && lhsType.getField(fieldName) != null)
            .toList();
        if( !overwritten.isEmpty() )
            addSoftError("Join fields " + overwritten + " are present in both records -- the left-hand values will be overwritten by the right-hand values", node);
    }

    /**
     * Get the fields that are returned directly by a key function,
     * e.g. `{ r -> r.id }` or `{ r -> tuple(r.id, r.type) }`.
     *
     * @param closure
     */
    private static Set<String> keyFunctionFields(ClosureExpression closure) {
        var param = closure.isParameterSpecified() ? closure.getParameters()[0].getName() : "it";
        var result = new HashSet<String>();
        var key = closureResult(closure);
        if( key == null )
            return result;
        var elements = key instanceof MethodCallExpression mce && "tuple".equals(mce.getMethodAsString())
            ? asMethodCallArguments(mce)
            : List.of(key);
        for( var element : elements ) {
            var field = fieldAccess(element, param);
            if( field != null )
                result.add(field);
        }
        return result;
    }

    /**
     * Get the fields that are compared directly by a join predicate,
     * e.g. `{ l, r -> l.id == r.id && l.type == r.type }`.
     *
     * @param closure
     */
    private static Set<String> joinPredicateFields(ClosureExpression closure) {
        var params = closure.getParameters();
        var result = new HashSet<String>();
        collectEqualFields(closureResult(closure), params[0].getName(), params[1].getName(), result);
        return result;
    }

    private static void collectEqualFields(Expression node, String lhsParam, String rhsParam, Set<String> result) {
        if( !(node instanceof BinaryExpression be) )
            return;
        var op = be.getOperation().getText();
        if( "&&".equals(op) ) {
            collectEqualFields(be.getLeftExpression(), lhsParam, rhsParam, result);
            collectEqualFields(be.getRightExpression(), lhsParam, rhsParam, result);
        }
        else if( "==".equals(op) ) {
            var a = be.getLeftExpression();
            var b = be.getRightExpression();
            var field = fieldAccess(a, lhsParam);
            if( field != null && field.equals(fieldAccess(b, rhsParam)) )
                result.add(field);
            field = fieldAccess(a, rhsParam);
            if( field != null && field.equals(fieldAccess(b, lhsParam)) )
                result.add(field);
        }
    }

    private static String fieldAccess(Expression node, String param) {
        if( node instanceof PropertyExpression pe && pe.getObjectExpression() instanceof VariableExpression ve && param.equals(ve.getName()) )
            return pe.getPropertyAsString();
        return null;
    }

    private static Expression closureResult(ClosureExpression closure) {
        var code = closure.getCode();
        var stmt = code instanceof BlockStatement block && !block.getStatements().isEmpty()
            ? block.getStatements().get(block.getStatements().size() - 1)
            : code;
        if( stmt instanceof ExpressionStatement es )
            return es.getExpression();
        if( stmt instanceof ReturnStatement rs )
            return rs.getExpression();
        return null;
    }

    private static ClassNode dataflowElementType(ClassNode type) {
        if( CHANNEL_TYPE.equals(type) || VALUE_TYPE.equals(type) )
            return elementType(type);
        return ClassHelper.dynamicType();
    }

    private void addError(String message, ASTNode node) {
        addError(message, node, false);
    }

    private void addSoftError(String message, ASTNode node) {
        addError(message, node, true);
    }

    private void addError(String message, ASTNode node, boolean softError) {
        var cause = new TypeError(message, node, softError);
        var errorMessage = new SyntaxErrorMessage(cause, sourceUnit);
        errorCollector.addErrorAndContinue(errorMessage);
    }

}
