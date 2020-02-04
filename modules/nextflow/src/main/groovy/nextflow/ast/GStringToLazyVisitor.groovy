/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.ast

import groovy.transform.CompileStatic
import org.codehaus.groovy.ast.ClassCodeVisitorSupport
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.GStringExpression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.control.SourceUnit

/**
 * Transform any GString to a Lazy GString i.e.
 *
 * from
 *   "${foo} ${bar}"
 * to
 *   "${->foo} ${->bar}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class GStringToLazyVisitor extends ClassCodeVisitorSupport {

    final SourceUnit sourceUnit

    private boolean withinClosure

    private List<String> names = []

    GStringToLazyVisitor(SourceUnit unit) {
        this.sourceUnit = unit
    }

    @Override
    void visitClosureExpression(ClosureExpression expression) {
        withinClosure = true
        try {
            super.visitClosureExpression(expression)
        }
        finally {
            withinClosure = false
        }
    }

    @Override
    void visitMethodCallExpression(MethodCallExpression call) {
        call.getObjectExpression().visit(this)
        call.getMethod().visit(this)
        names.push(call.methodAsString)
        call.getArguments().visit(this)
        names.pop()
    }

    @Override
    void visitGStringExpression(GStringExpression expression) {
        // channels values are not supposed to be lazy evaluated
        // therefore stop visiting when reaching `from` keyword
        if( !withinClosure && names.last()!='from' ) {
            xformToLazy(expression)
        }
    }

    protected void xformToLazy(GStringExpression str) {
        def values = str.getValues()
        def normalised = new Expression[values.size()]

        // wrap all non-closure to a ClosureExpression
        for( int i=0; i<values.size(); i++ ) {
            final item = values[i]
            if( item instanceof ClosureExpression  ) {
                // when there is already a closure the conversion it is skipped
                // because it's supposed the gstring is already a lazy-string
                return
            }
            normalised[i] = wrapWithClosure(item)
        }

        for( int i=0; i<values.size(); i++ ) {
            values[i] = normalised[i]
        }
    }

    protected ClosureExpression wrapWithClosure( Expression expr ) {

        // create an expression statement for the given `expression`
        def statement = new ExpressionStatement(expr)
        // add it to a new block
        def block = new BlockStatement()
        block.addStatement(statement)
        // create a closure over the given block
        // note: the closure parameter argument must be *null* to force the creation of a closure like {-> something}
        // otherwise it creates a closure with an implicit parameter that is managed in a different manner by the
        // GString -- see http://docs.groovy-lang.org/latest/html/documentation/#_special_case_of_interpolating_closure_expressions
        new ClosureExpression( null, block )

    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit
    }

}
