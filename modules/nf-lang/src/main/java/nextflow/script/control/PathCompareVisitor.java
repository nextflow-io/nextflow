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

package nextflow.script.control;

import org.codehaus.groovy.ast.ClassCodeExpressionTransformer;
import org.codehaus.groovy.ast.expr.BinaryExpression;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.syntax.Types;

import static nextflow.script.ast.ASTUtils.*;
import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Replace path comparisons with explicit method calls.
 *
 * This is required to correctly compare `Path` values, which are
 * not supported by default because `Path` implements `Comparable`.
 *
 * @see https://stackoverflow.com/questions/28355773/in-groovy-why-does-the-behaviour-of-change-for-interfaces-extending-compar#comment45123447_28387391
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class PathCompareVisitor extends ClassCodeExpressionTransformer {

    private SourceUnit sourceUnit;

    public PathCompareVisitor(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    @Override
    public Expression transform(Expression node) {
        if( node instanceof BinaryExpression be ) {
            return transformBinaryExpression(be);
        }

        if( node instanceof ClosureExpression ce ) {
            ce.visit(this);
            return ce;
        }

        return super.transform(node);
    }

    /**
     * Replace comparison operations with explicit calls to
     * the appropriate {@link LangHelpers} method.
     *
     * @param node
     */
    protected Expression transformBinaryExpression(BinaryExpression node) {

        var left = node.getLeftExpression();
        var right = node.getRightExpression();

        return switch( node.getOperation().getType() ) {
            case Types.COMPARE_EQUAL ->
                call("compareEqual", left, right);

            case Types.COMPARE_NOT_EQUAL ->
                notX(call("compareEqual", left, right));

            case Types.COMPARE_LESS_THAN ->
                call("compareLessThan", left, right);

            case Types.COMPARE_LESS_THAN_EQUAL ->
                call("compareLessThanEqual", left, right);

            case Types.COMPARE_GREATER_THAN ->
                call("compareGreaterThan", left, right);

            case Types.COMPARE_GREATER_THAN_EQUAL ->
                call("compareGreaterThanEqual", left, right);

            default -> super.transform(node);
        };
    }

    private MethodCallExpression call(String method, Expression left, Expression right) {
        return callX(
            classX("nextflow.util.LangHelpers"),
            method,
            args(transform(left), transform(right)));
    }

}
