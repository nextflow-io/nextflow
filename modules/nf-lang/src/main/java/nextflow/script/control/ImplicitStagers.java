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
package nextflow.script.control;

import java.nio.file.Path;
import java.util.Collection;

import nextflow.script.ast.RecordNode;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.Variable;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;

import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Inference of the implicit staging directives that a typed input declaration
 * carries: which declarations mean "stage this into the task directory".
 *
 * Shared by {@link ProcessToGroovyVisitorV2} and {@link AgentToGroovyVisitor} so
 * that the same declaration means the same thing in a process and in an agent —
 * the element type of a {@code Collection<Path>} and the fields of a record type
 * are only visible at compile time, so this cannot be re-derived at runtime.
 */
class ImplicitStagers {

    /**
     * Add the implicit staging directives inferred from the declared input type:
     *
     * - Inputs with type Path or a Path collection (e.g. Set&lt;Path&gt;)
     *   are staged as input files.
     *
     * - Inputs with a record type are recursively inspected for nested
     *   file inputs based on the record type definition.
     *
     * - Inputs with a record collection type (e.g. List&lt;Sample&gt;) are
     *   inspected in the same way, staging the files of every element.
     *
     * @param param
     * @param target
     * @param stagers
     */
    static void visitInputType(Variable param, Expression target, BlockStatement stagers) {
        visitInputType(param.getType(), target, false, stagers);
    }

    private static void visitInputType(ClassNode cn, Expression target, boolean spread, BlockStatement stagers) {
        if( isPathType(cn) ) {
            var value = spread ? flattenNonNull(target) : target;
            var stager = stmt(callThisX("stageAs", args(closureX(stmt(value)))));
            stagers.addStatement(stager);
        }
        else if( isRecordType(cn) ) {
            for( var fn : cn.getFields() ) {
                var property = propX(target, fn.getName());
                property.setSpreadSafe(spread);
                visitInputType(fn.getType(), property, spread, stagers);
            }
        }
        else if( isRecordCollectionType(cn) ) {
            var elements = spread ? safeCallX(target, "flatten") : target;
            visitInputType(cn.getGenericsTypes()[0].getType(), elements, true, stagers);
        }
    }

    /**
     * Flatten the values collected from a record collection and
     * discard null values (i.e. from nullable fields).
     *
     * @param target
     */
    private static Expression flattenNonNull(Expression target) {
        var p = param(ClassHelper.dynamicType(), "p");
        var notNull = closureX(params(p), stmt(notNullX(varX(p))));
        return safeCallX(safeCallX(target, "flatten"), "findAll", args(notNull));
    }

    private static Expression safeCallX(Expression target, String name) {
        return safeCallX(target, name, MethodCallExpression.NO_ARGUMENTS);
    }

    private static Expression safeCallX(Expression target, String name, Expression args) {
        var call = callX(target, name, args);
        call.setSafe(true);
        return call;
    }

    private static boolean isPathType(ClassNode cn) {
        if( !cn.isResolved() )
            return false;
        var clazz = cn.getTypeClass();
        if( Path.class.isAssignableFrom(clazz) ) {
            return true;
        }
        if( Collection.class.isAssignableFrom(clazz) && cn.isUsingGenerics() ) {
            var elementType = cn.getGenericsTypes()[0].getType();
            return Path.class.isAssignableFrom(elementType.getTypeClass());
        }
        return false;
    }

    private static boolean isRecordType(ClassNode cn) {
        return cn.redirect() instanceof RecordNode;
    }

    private static boolean isRecordCollectionType(ClassNode cn) {
        if( !cn.isResolved() || !cn.isUsingGenerics() )
            return false;
        if( !Collection.class.isAssignableFrom(cn.getTypeClass()) )
            return false;
        return isRecordType(cn.getGenericsTypes()[0].getType());
    }

}
