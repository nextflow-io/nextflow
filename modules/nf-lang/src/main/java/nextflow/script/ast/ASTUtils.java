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

package nextflow.script.ast;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Stream;

import groovy.transform.NamedParams;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.AnnotatedNode;
import org.codehaus.groovy.ast.AnnotationNode;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.PropertyNode;
import org.codehaus.groovy.ast.Variable;
import org.codehaus.groovy.ast.expr.AnnotationConstantExpression;
import org.codehaus.groovy.ast.expr.ClassExpression;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.ConstantExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.ListExpression;
import org.codehaus.groovy.ast.expr.MapExpression;
import org.codehaus.groovy.ast.expr.MapEntryExpression;
import org.codehaus.groovy.ast.expr.MethodCall;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.NamedArgumentListExpression;
import org.codehaus.groovy.ast.expr.PropertyExpression;
import org.codehaus.groovy.ast.expr.TupleExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.ast.tools.GeneralUtils;

/**
 * Utility functions for common AST operations.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ASTUtils {

    public static Expression classX(String name) {
        return GeneralUtils.classX(ClassHelper.makeWithoutCaching(name));
    }

    public static Expression createX(Class type, TupleExpression args) {
        return GeneralUtils.ctorX(new ClassNode(type), args);
    }

    public static Expression createX(Class type, Expression... expressions) {
        return GeneralUtils.ctorX(new ClassNode(type), GeneralUtils.args(expressions));
    }

    public static Expression createX(String name, TupleExpression args) {
        return GeneralUtils.ctorX(ClassHelper.makeWithoutCaching(name), args);
    }

    public static Expression createX(String name, Expression... expressions) {
        return GeneralUtils.ctorX(ClassHelper.makeWithoutCaching(name), GeneralUtils.args(expressions));
    }

    public static List<Statement> asBlockStatements(Statement statement) {
        return statement instanceof BlockStatement block
            ? block.getStatements()
            : Collections.emptyList();
    }

    public static ConstantExpression asConstX(Expression expression) {
        return expression instanceof ConstantExpression ce ? ce : null;
    }

    /**
     * Given a statement which represents a block statement of directives,
     * iterate through each directive as a method call expression.
     *
     * @param statement
     */
    public static Stream<MethodCallExpression> asDirectives(Statement statement) {
        return asBlockStatements(statement)
            .stream()
            .map(stmt -> asMethodCallX(stmt))
            .filter(mce -> mce != null);
    }

    /**
     * Given a method call which represents a definition (i.e. DSL) block, get
     * the definition body, which is the block statement of the last closure argument
     * in the method call.
     *
     * @param call
     * @param argsCount
     */
    public static BlockStatement asDslBlock(MethodCallExpression call, int argsCount) {
        var args = asMethodCallArguments(call);
        if( args.size() != argsCount )
            return null;
        var lastArg = args.get(args.size() - 1);
        if( !(lastArg instanceof ClosureExpression) )
            return null;
        var closure = (ClosureExpression) lastArg;
        return (BlockStatement) closure.getCode();
    }

    public static Parameter[] asFlatParams(Parameter[] params) {
        return Arrays.stream(params)
            .flatMap((param) -> (
                param instanceof TupleParameter tp
                    ? Arrays.stream(tp.components)
                    : Stream.of(param)
            ))
            .toArray(Parameter[]::new);
    }

    public static MethodCallExpression asMethodCallX(Statement stmt) {
        if( !(stmt instanceof ExpressionStatement) )
            return null;
        var stmtX = (ExpressionStatement) stmt;
        if( !(stmtX.getExpression() instanceof MethodCallExpression) )
            return null;
        return (MethodCallExpression) stmtX.getExpression();
    }

    public static List<Expression> asMethodCallArguments(MethodCall call) {
        return ((TupleExpression) call.getArguments()).getExpressions();
    }

    public static List<MapEntryExpression> asNamedArgs(MethodCall call) {
        var args = asMethodCallArguments(call);
        return args.size() > 0 && args.get(0) instanceof NamedArgumentListExpression nale
            ? nale.getMapEntryExpressions()
            : Collections.emptyList();
    }

    /**
     * Given a parameter with a @NamedParams annotation,
     * return the map of named params.
     *
     * @param parameter
     */
    public static Map<String, AnnotationNode> asNamedParams(Parameter parameter) {
        var namedParams = new LinkedHashMap<String, AnnotationNode>();
        parameter.getAnnotations().stream()
            .filter(an -> an.getClassNode().getName().equals(NamedParams.class.getName()))
            .flatMap(an -> {
                var value = an.getMember("value");
                return value instanceof ListExpression le
                    ? le.getExpressions().stream()
                    : Stream.empty();
            })
            .forEach((value) -> {
                if( !(value instanceof AnnotationConstantExpression) )
                    return;
                var ace = (AnnotationConstantExpression) value;
                var namedParam = (AnnotationNode) ace.getValue();
                var name = namedParam.getMember("value").getText();
                namedParams.put(name, namedParam);
            });
        return namedParams;
    }

    public static Parameter asNamedParam(AnnotationNode node) {
        var name = node.getMember("value").getText();
        var typeX = (ClassExpression) node.getMember("type");
        var type = typeX != null ? typeX.getType() : ClassHelper.dynamicType();
        return new Parameter(type, name);
    }

    public static VariableExpression asVarX(Statement statement) {
        return statement instanceof ExpressionStatement es ? asVarX(es.getExpression()) : null;
    }

    public static VariableExpression asVarX(Expression expression) {
        return expression instanceof VariableExpression ve ? ve : null;
    }

    /**
     * Given a variable which represents a method being accessed
     * as a variable, return the underlying method.
     *
     * @param variable
     */
    public static MethodNode asMethodVariable(Variable variable) {
        if( variable instanceof PropertyNode pn ) {
            if( pn.getNodeMetaData(ASTNodeMarker.METHOD_VARIABLE_TARGET) instanceof MethodNode mn )
                return mn;
        }
        return null;
    }

    /**
     * Given a property expression which represents a process or workflow
     * output, return the underlying process or workflow.
     *
     * @param node
     */
    public static MethodNode asMethodOutput(PropertyExpression node) {
        if( node.getObjectExpression() instanceof VariableExpression ve && "out".equals(node.getPropertyAsString()) )
            return asMethodVariable(ve.getAccessedVariable());
        return null;
    }

    /**
     * Given an annotated node (e.g. class, field, method), Find the first
     * annotation of the given type in the node's list of annotations.
     *
     * @param node
     * @param type
     */
    public static Optional<AnnotationNode> findAnnotation(AnnotatedNode node, Class type) {
        return node.getAnnotations().stream()
            .filter(an -> an.getClassNode().getName().equals(type.getName()))
            .findFirst();
    }
}
