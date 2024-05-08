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
 *
 */

package nextflow.config

import static org.codehaus.groovy.ast.tools.GeneralUtils.*

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassCodeExpressionTransformer
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.expr.ArgumentListExpression
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.ConstantExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.PropertyExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.control.CompilePhase
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.transform.ASTTransformation
import org.codehaus.groovy.transform.GroovyASTTransformation
/**
 * AST transformation that replaces properties prefixed with `secrets.`
 * with a static string
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@GroovyASTTransformation(phase = CompilePhase.CONVERSION)
class StripSecretsXformImpl implements ASTTransformation {

    private SourceUnit unit

    @Override
    void visit(ASTNode[] nodes, SourceUnit source) {
        this.unit = source
        createVisitor().visitClass((ClassNode)nodes[1])
    }

    protected ClassCodeExpressionTransformer createVisitor() {

        new ClassCodeExpressionTransformer() {

            private boolean isIncludeConfigArgument

            protected SourceUnit getSourceUnit() { unit }

            @Override
            Expression transform(Expression expr) {
                if (expr == null)
                    return null

                final result = replaceProperty(expr)
                if( result ) {
                    return result
                }
                if( expr instanceof MethodCallExpression ) {
                    return transformMethodCall(expr as MethodCallExpression)
                }
                if( expr instanceof ClosureExpression) {
                    visitClosureExpression(expr)
                }
                return super.transform(expr)
            }

            protected Expression transformMethodCall(MethodCallExpression call) {
                if( call.methodAsString=='includeConfig' && call.arguments instanceof ArgumentListExpression ) {
                    // flag method calls to `includeConfig` to preserve the use of secret properties
                    // for context see https://github.com/nextflow-io/nextflow/pull/4177
                    isIncludeConfigArgument = true
                    try {
                        return super.transform(call)
                    }
                    finally {
                        isIncludeConfigArgument = false
                    }
                }
                else {
                    super.transform(call)
                }
            }

            protected Expression replaceProperty(Expression expr) {
                if( expr !instanceof PropertyExpression )
                    return null
                final p = expr as PropertyExpression
                final isSecretProperty =
                    p.objectExpression instanceof VariableExpression
                        && p.objectExpression.text=='secrets'
                        && p.property instanceof ConstantExpression
                if( isSecretProperty && !isIncludeConfigArgument ) {
                    // replace the reference to a secrets property into a constant
                    // in order to disclose the secret value when using the `nextflow config` command
                    return constX("secrets.${p.propertyAsString}".toString())
                }
                return null
            }

        }
    }

}
