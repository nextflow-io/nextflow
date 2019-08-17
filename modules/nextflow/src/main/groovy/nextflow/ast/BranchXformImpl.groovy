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
import groovy.util.logging.Slf4j
import nextflow.script.TokenBranchChoice
import nextflow.script.TokenBranchDef
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassCodeExpressionTransformer
import org.codehaus.groovy.ast.ClassCodeVisitorSupport
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.Parameter
import org.codehaus.groovy.ast.VariableScope
import org.codehaus.groovy.ast.expr.ArgumentListExpression
import org.codehaus.groovy.ast.expr.BooleanExpression
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.EmptyStatement
import org.codehaus.groovy.ast.stmt.IfStatement
import org.codehaus.groovy.ast.stmt.ReturnStatement
import org.codehaus.groovy.ast.stmt.Statement
import org.codehaus.groovy.ast.tools.GeneralUtils
import org.codehaus.groovy.control.CompilePhase
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.transform.ASTTransformation
import org.codehaus.groovy.transform.GroovyASTTransformation
import static nextflow.ast.ASTHelpers.createX
import static nextflow.ast.ASTHelpers.isArgsX
import static nextflow.ast.ASTHelpers.isBlockStmt
import static nextflow.ast.ASTHelpers.isClosureX
import static nextflow.ast.ASTHelpers.isStmtX
import static nextflow.ast.ASTHelpers.syntaxError
import static org.codehaus.groovy.ast.tools.GeneralUtils.block
import static org.codehaus.groovy.ast.tools.GeneralUtils.closureX
import static org.codehaus.groovy.ast.tools.GeneralUtils.constX
import static org.codehaus.groovy.ast.tools.GeneralUtils.listX
import static org.codehaus.groovy.ast.tools.GeneralUtils.returnS
import static org.codehaus.groovy.ast.tools.GeneralUtils.stmt
import static org.codehaus.groovy.ast.tools.GeneralUtils.varX
/**
 * Implements Nextflow operator xform logic
 * See http://groovy-lang.org/metaprogramming.html#_classcodeexpressiontransformer
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@GroovyASTTransformation(phase = CompilePhase.CONVERSION)
class BranchXformImpl implements ASTTransformation {

    static final public String OPERATOR_NAME = 'branch'

    static final public String CRITERIA_NAME = 'branchCriteria'

    SourceUnit unit

    static class BranchCondition {
        String label
        Expression condition
        List<Statement> code = new ArrayList<>()

        BranchCondition(String lbl, Expression condition) {
            this.label = lbl
            this.condition = condition
        }
    }

    @CompileStatic
    class Transformer {

        final List<Statement> impl = new ArrayList<>(20)
        final Set<String> allLabels = new LinkedHashSet<>(10)
        final Map<String,BranchCondition> allBlocks = new LinkedHashMap<>(10)
        final ClosureExpression body
        final BlockStatement code
        final VariableScope scope
        final MethodCallExpression method
        String current = null

        Transformer(MethodCallExpression method, ClosureExpression body) {
            this.body = body
            this.code = isBlockStmt(body.code)
            this.scope = body.variableScope
            this.method = method
        }

        protected Expression apply() {
            if( !code )
                return method

            if( body.parameters==null )
                syntaxError("Branch evaluation closure should declare at least one parameter or use the implicit `it` parameter", body, unit)

            for( Statement stmt : code.statements ) {
                final label = stmt.statementLabel
                if( label ) {
                    if( !allLabels.add(label) ) {
                        syntaxError("Branch label already used: $label: $stmt", stmt, unit)
                        break
                    }
                    current = label

                    final expr = isStmtX(stmt)?.expression
                    final block = new BranchCondition(label, expr)
                    allBlocks.put(label, block)

                    if(!block.condition)
                        syntaxError("Unexpected statement here: $label: $stmt", stmt, unit)

                }
                else if( current ) {
                    allBlocks.get(current).code.add(stmt)
                }
                else {
                    impl.add(stmt)
                }
            }

            if( !allBlocks ) {
                syntaxError("Branch evaluation closure should contain at least one branch expression", body, unit)
                return null
            }

            // now write the code
            IfStatement previousIf
            IfStatement rootIf = null
            IfStatement currentIf = null
            for( BranchCondition entry : allBlocks.values()) {
                previousIf = currentIf
                // the if block
                final cond = new BooleanExpression(entry.condition)
                currentIf = new IfStatement(cond, blockS(entry), null)
                currentIf.statementLabel = entry.label
                if( previousIf )
                    previousIf.elseBlock = currentIf
                if( !rootIf )
                    rootIf = currentIf
            }

            if( !currentIf.elseBlock )
                currentIf.elseBlock = new EmptyStatement()

            // add the main if statement to the implementation statements list
            impl.add(rootIf)

            // finally create a new closure with all the implementation code
            def main = closureX(body.parameters, block(scope,impl))
            def result = createX(TokenBranchDef, main, GeneralUtils.list2args(new ArrayList(allLabels)))

            // create the implementation closure
            final closure = new ClosureExpression(null, block(scope, stmt(result)))
            // pass the new closure as the method argument
            method.setArguments( new ArgumentListExpression(closure) )
            // return the modified method call
            return method
        }

        protected Expression paramX(Parameter[] params) {
            assert params!=null
            if( params.size()==0 ) {
               return varX('it')
            }

            if( params.size()==1 ) {
                return varX(params[0].name)
            }

            final ret = new ArrayList<Expression>(params.size())
            for( Parameter p : params ) {
                ret.add( varX(p.name) )
            }
            return listX(ret)
        }

        protected BlockStatement blockS(BranchCondition branch) {
            List<Statement> statements = branch.code
            if( !statements ) {
                // no result --> return the implicit closure param
                statements.add( resultS(paramX(body.parameters), branch.label) )
                return block(scope,statements)
            }

            // otherwise check for the first return statement and
            if( fixReturnStatement(statements, branch.label)) {
                return block(scope,statements)
            }

            // if there's no return statement add it using the last expression
            final last = statements.size()-1
            final stmt = isStmtX(statements[last])
            if( stmt ) {
                statements[last] = resultS(stmt.expression, branch.label)
                return block(scope,statements)
            }

            syntaxError("Unexpected statement in branch condition", statements[last], unit)
            return null
        }

        protected Statement resultS(Expression expr, String choice) {
            returnS(
                    createX( TokenBranchChoice, expr, constX(choice)
            ))
        }

        protected boolean fixReturnStatement(List<Statement> statements, String choice) {

            final resultXform = new ClassCodeVisitorSupport() {

                int count

                @Override
                protected SourceUnit getSourceUnit() { unit }

                @Override
                void visitReturnStatement(ReturnStatement statement) {
                    count++
                    statement.expression = createX(TokenBranchChoice, statement.expression, constX(choice))
                    super.visitReturnStatement(statement)
                }

                int apply() {
                    for (Statement statement : statements) {
                        statement.visit(this);
                    }

                    return count
                }
            }

            resultXform.apply()>0
        }
    }

    @Override
    void visit(ASTNode[] nodes, SourceUnit source) {
        this.unit = unit
        createVisitor().visitClass((ClassNode)nodes[1])
    }

    protected ClosureExpression isSwitchOpCall(Expression expr) {
        final m = ASTHelpers.isMethodCallX(expr)
        if( m ) {
            final name = m.methodAsString
            final args = isArgsX(m.arguments)
            final ClosureExpression clo = args && args.size()>0 ? isClosureX(args.last()) : null
            return ((name==OPERATOR_NAME && args.size()==1) ||
                    (name==CRITERIA_NAME && args.size()==1 ) ? clo : null)
        }
        return null
    }


    protected ClassCodeExpressionTransformer createVisitor() {

        new ClassCodeExpressionTransformer() {

            protected SourceUnit getSourceUnit() { unit }

            @Override
            Expression transform(Expression expr) {
                if (expr == null)
                    return null

                final body = isSwitchOpCall(expr)
                if( body ) {
                    return new Transformer(expr as MethodCallExpression, body).apply()
                }
                else if( expr instanceof ClosureExpression) {
                    visitClosureExpression(expr)
                }

                return super.transform(expr)
            }
        }
    }

}


