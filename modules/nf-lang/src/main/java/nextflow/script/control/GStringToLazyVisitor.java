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

import org.codehaus.groovy.ast.ClassCodeVisitorSupport;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.GStringExpression;
import org.codehaus.groovy.control.SourceUnit;

import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Transform a GString to a Lazy GString.
 *
 * from
 *   "${foo} ${bar}"
 * to
 *   "${->foo} ${->bar}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class GStringToLazyVisitor extends ClassCodeVisitorSupport {

    private SourceUnit sourceUnit;

    private boolean inClosure;

    public GStringToLazyVisitor(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    @Override
    public void visitClosureExpression(ClosureExpression node) {
        inClosure = true;
        try {
            super.visitClosureExpression(node);
        }
        finally {
            inClosure = false;
        }
    }

    @Override
    public void visitGStringExpression(GStringExpression node) {
        // gstrings in a closure will be lazily evaluated and therefore
        // don't need to be lazy themselves
        if( !inClosure ) {
            transformToLazy(node);
        }
    }

    private void transformToLazy(GStringExpression node) {
        var values = node.getValues();
        var lazyValues = new Expression[values.size()];

        // wrap all non-closure expressions in a closure
        for( int i = 0; i < values.size(); i++ ) {
            var value = values.get(i);
            if( value instanceof ClosureExpression ) {
                // when the value is already a closure, skip the entire gstring
                // because it is assumed to already be lazy
                return;
            }
            lazyValues[i] = wrapExpressionInClosure(value);
        }

        for( int i = 0; i < values.size(); i++ ) {
            values.set(i, lazyValues[i]);
        }
    }

    protected ClosureExpression wrapExpressionInClosure(Expression node)  {
        // note: the closure parameter argument must be *null* to force the creation of a closure like {-> something}
        // otherwise it creates a closure with an implicit parameter that is managed in a different manner by the
        // GString -- see http://docs.groovy-lang.org/latest/html/documentation/#_special_case_of_interpolating_closure_expressions
        return closureX(null, block(stmt(node)));
    }

}
