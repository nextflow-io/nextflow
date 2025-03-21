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
package nextflow.script.types;

import java.util.Collections;
import java.util.List;

import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.ast.ProcessNode;
import nextflow.script.ast.WorkflowNode;
import nextflow.script.dsl.Constant;
import nextflow.script.dsl.Description;
import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.FieldNode;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.Variable;
import org.codehaus.groovy.ast.expr.ClassExpression;
import org.codehaus.groovy.ast.expr.ConstructorCallExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.MethodCall;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.PropertyExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;

import static nextflow.script.ast.ASTUtils.*;

/**
 * Utility methods for inferring and checking the type of an expression.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class TypeChecker {

    /**
     * Get the type (i.e. class node) of a variable or expression.
     *
     * @param node
     */
    public static ClassNode getType(ASTNode node) {
        if( node.getNodeMetaData(ASTNodeMarker.INFERRED_TYPE) instanceof ClassNode cn )
            return cn;
        var inferredType = inferredType(node);
        var result = Types.SHIM_TYPES.containsKey(inferredType)
            ? Types.SHIM_TYPES.get(inferredType)
            : inferredType;
        node.putNodeMetaData(ASTNodeMarker.INFERRED_TYPE, result);
        return result;
    }

    private static ClassNode inferredType(ASTNode node) {
        if( node instanceof ClassExpression ce ) {
            return ce.getType();
        }

        if( node instanceof ConstructorCallExpression cce ) {
            return cce.getType();
        }

        if( node instanceof MethodCallExpression mce ) {
            var mn = inferMethodTarget(mce);
            return mn != null ? mn.getReturnType() : null;
        }

        if( node instanceof PropertyExpression pe ) {
            var mn = asMethodOutput(pe);
            if( mn instanceof ProcessNode pn )
                return pn.getReturnType();
            if( mn instanceof WorkflowNode wn )
                return wn.getReturnType();
            var fn = inferPropertyTarget(pe);
            return fn != null ? getType(fn) : null;
        }

        if( node instanceof VariableExpression ve ) {
            var mn = asMethodVariable(ve.getAccessedVariable());
            if( mn instanceof ProcessNode || mn instanceof WorkflowNode )
                return null;
        }

        if( node instanceof Variable variable && variable.getOriginType() != null ) {
            return variable.getOriginType();
        }

        return node instanceof Expression exp
            ? exp.getType()
            : null;
    }

    /**
     * Find the field node targeted by a property expression.
     *
     * @param node
     */
    public static FieldNode inferPropertyTarget(PropertyExpression node) {
        var receiverType = getType(node.getObjectExpression());
        if( receiverType == null )
            return null;
        var name = node.getPropertyAsString();
        var result = receiverType.getField(name);
        if( result != null )
            return result;
        return receiverType.getMethods().stream()
            .filter((mn) -> {
                if( !mn.isPublic() )
                    return false;
                var an = findAnnotation(mn, Constant.class);
                if( !an.isPresent() )
                    return false;
                return name.equals(an.get().getMember("value").getText());
            })
            .map((mn) -> {
                var fn = new FieldNode(name, 0xF, mn.getReturnType(), mn.getDeclaringClass(), null);
                findAnnotation(mn, Description.class).ifPresent((an) -> {
                    fn.addAnnotation(an);
                });
                return fn;
            })
            .findFirst().orElse(null);
    }

    /**
     * Find the method node targeted by a call expression.
     *
     * @param node
     */
    public static MethodNode inferMethodTarget(MethodCall node) {
        var methods = methodOverloads(node);
        var arguments = asMethodCallArguments(node);
        return methods.stream()
            .filter((mn) -> {
                var parameters = mn.getParameters();
                if( arguments.size() != parameters.length ) {
                    return false;
                }
                for( int i = 0; i < parameters.length; i++ ) {
                    var paramType = parameters[i].getType();
                    var argType = arguments.get(i).getType();
                    if( !Types.isAssignableFrom(paramType, argType) )
                        return false;
                }
                return true;
            })
            .findFirst()
            .orElse(null);
    }

    private static List<MethodNode> methodOverloads(MethodCall node) {
        if( node instanceof MethodCallExpression mce ) {
            if( mce.getNodeMetaData(ASTNodeMarker.METHOD_TARGET) instanceof MethodNode mn )
                return List.of(mn);

            if( !mce.isImplicitThis() ) {
                var receiverType = getType(mce.getObjectExpression());
                if( receiverType != null )
                    return methodsForType(receiverType, mce.getMethodAsString());
            }
        }

        if( node instanceof ConstructorCallExpression cce ) {
            var constructorType = cce.getType();
            if( constructorType != null ) {
                return constructorType.getDeclaredConstructors().stream()
                    .map(ctor -> (MethodNode) ctor)
                    .toList();
            }
        }

        return Collections.emptyList();
    }

    private static List<MethodNode> methodsForType(ClassNode cn, String name) {
        try {
            return cn.getAllDeclaredMethods().stream()
                .filter(mn -> mn.getName().equals(name))
                .filter(mn -> !ClassHelper.isObjectType(mn.getDeclaringClass()))
                .toList();
        }
        catch( NullPointerException e ) {
            return Collections.emptyList();
        }
    }

}
