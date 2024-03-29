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

package nextflow.ast

import groovy.transform.CompileStatic
import org.codehaus.groovy.ast.ClassCodeVisitorSupport
import org.codehaus.groovy.ast.expr.ArgumentListExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.MapExpression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.ast.stmt.Statement
import org.codehaus.groovy.control.SourceUnit

import static org.codehaus.groovy.ast.tools.GeneralUtils.*

/**
 * Extract unstaging directives from process output
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessOutputVisitor extends ClassCodeVisitorSupport {

    final SourceUnit sourceUnit

    final List<Statement> statements = []

    private int evalCount = 0

    private int pathCount = 0

    ProcessOutputVisitor(SourceUnit unit) {
        this.sourceUnit = unit
    }

    @Override
    void visitMethodCallExpression(MethodCallExpression call) {
        extractDirective(call)
        super.visitMethodCallExpression(call)
    }

    void extractDirective(MethodCallExpression call) {
        if( call.objectExpression?.text != 'this' )
            return

        if( call.arguments !instanceof ArgumentListExpression )
            return

        final name = call.methodAsString
        final args = (ArgumentListExpression)call.arguments

        /**
         * env(name) -> _typed_out_env(name)
         */
        if( name == 'env' ) {
            if( args.size() != 1 )
                return

            statements << makeDirective(
                '_typed_out_env',
                args[0]
            )
        }

        /**
         * eval(cmd) -> _typed_out_eval(key) { cmd }
         *           -> eval(key)
         */
        else if( name == 'eval' ) {
            if( args.size() != 1 )
                return

            final key = constX("nxf_out_eval_${evalCount++}".toString())

            statements << makeDirective(
                '_typed_out_eval',
                key,
                closureX(new ExpressionStatement(args[0]))
            )

            call.arguments = new ArgumentListExpression(key)
        }

        /**
         * path(opts, pattern) -> _typed_out_path(opts, key) { pattern }
         *                     -> path(key)
         */
        else if( name == 'path' ) {
            def opts = null
            def pattern
            if( args.size() == 1 ) {
                pattern = args[0]
            }
            else if( args.size() == 2 ) {
                opts = args[0]
                pattern = args[1]
            }
            else return

            final key = constX("\$file${pathCount++}".toString())

            statements << makeDirective(
                '_typed_out_path',
                opts ?: new MapExpression(),
                key,
                closureX(new ExpressionStatement(pattern))
            )

            call.arguments = new ArgumentListExpression(key)
        }
    }

    Statement makeDirective(String name, Expression... args) {
        new ExpressionStatement(
            new MethodCallExpression(
                new VariableExpression('this'),
                name,
                new ArgumentListExpression(args)
            )
        )
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit
    }

}
