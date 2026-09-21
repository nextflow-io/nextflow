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
import java.util.List;

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
import org.codehaus.groovy.ast.expr.ConstantExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.NamedArgumentListExpression;
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
        checkJoinByFields(lhsType, rhsType, arguments);
        var elementType = recordSumType(lhsType, rhsType);
        return makeType(CHANNEL_TYPE, elementType);
    }

    private void checkJoinByFields(ClassNode lhsType, ClassNode rhsType, List<Expression> arguments) {
        // Report an error if the `by` field of a `join` operation is not present in both records.
        if( arguments.isEmpty() || !(arguments.get(0) instanceof NamedArgumentListExpression nale) )
            return;
        for( var entry : nale.getMapEntryExpressions() ) {
            if( !"by".equals(entry.getKeyExpression().getText()) )
                continue;
            var field = entry.getValueExpression();
            if( !(field instanceof ConstantExpression ce) || !(ce.getValue() instanceof String name) )
                continue;
            if( lhsType.getField(name) == null )
                addError("Join field `" + name + "` is not present in left-hand side", field);
            if( rhsType.getField(name) == null )
                addError("Join field `" + name + "` is not present in right-hand side", field);
        }
    }

    private static ClassNode dataflowElementType(ClassNode type) {
        if( CHANNEL_TYPE.equals(type) || VALUE_TYPE.equals(type) )
            return elementType(type);
        return ClassHelper.dynamicType();
    }

    private void addError(String message, ASTNode node) {
        var cause = new TypeError(message, node);
        var errorMessage = new SyntaxErrorMessage(cause, sourceUnit);
        errorCollector.addErrorAndContinue(errorMessage);
    }

}
