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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.ClassCodeExpressionTransformer;
import org.codehaus.groovy.ast.CodeVisitorSupport;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.TupleExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.EmptyStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.IfStatement;
import org.codehaus.groovy.ast.stmt.ReturnStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.runtime.DefaultGroovyMethods;
import org.codehaus.groovy.syntax.SyntaxException;

import static nextflow.script.ast.ASTUtils.*;
import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Transform closure arguments for branch and multiMap
 * into the appropriate criteria objects.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class OpCriteriaVisitor extends ClassCodeExpressionTransformer {

    private SourceUnit sourceUnit;

    public OpCriteriaVisitor(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    @Override
    public Expression transform(Expression node) {
        if( node instanceof ClosureExpression ce )
            super.visitClosureExpression(ce);

        if( node instanceof MethodCallExpression mce )
            return transformMethodCall(mce);

        return super.transform(node);
    }

    private Expression transformMethodCall(MethodCallExpression node) {
        ClosureExpression body;
        if( (body=asBranchOpClosure(node)) != null ) {
            var criteria = new BranchTransformer(body).apply();
            node.setArguments(args(criteria));
            return node;
        }

        if( (body=asMultiMapOpClosure(node)) != null ) {
            var criteria = new MultiMapTransformer(body).apply();
            node.setArguments(args(criteria));
            return node;
        }

        return super.transform(node);
    }

    private ClosureExpression asBranchOpClosure(MethodCallExpression node) {
        var name = node.getMethodAsString();
        var arguments = (TupleExpression) node.getArguments();
        var argCount = arguments.getExpressions().size();
        var closureArg = argCount > 0 && arguments.getExpression(argCount - 1) instanceof ClosureExpression ce ? ce : null;
        if( "branch".equals(name) && argCount == 1 )
            return closureArg;
        if( node.isImplicitThis() && "branchCriteria".equals(name) && argCount == 1 )
            return closureArg;
        return null;
    }

    private ClosureExpression asMultiMapOpClosure(MethodCallExpression node) {
        var name = node.getMethodAsString();
        var arguments = (TupleExpression) node.getArguments();
        var argCount = arguments.getExpressions().size();
        var closureArg = argCount > 0 && arguments.getExpression(argCount - 1) instanceof ClosureExpression ce ? ce : null;
        if( "multiMap".equals(name) && argCount == 1 )
            return closureArg;
        if( node.isImplicitThis() && "multiMapCriteria".equals(name) && argCount == 1 )
            return closureArg;
        return null;
    }

    private void syntaxError(String message, ASTNode node) {
        sourceUnit.addError(new SyntaxException(message, node.getLineNumber(), node.getColumnNumber()));
    }

    class BranchTransformer {

        private ClosureExpression body;

        public BranchTransformer(ClosureExpression body) {
            this.body = body;
        }

        Expression apply() {
            if( body.getParameters() == null ) {
                syntaxError("Branch criteria should declare at least one parameter or use the implicit `it` parameter", body);
                return body;
            }

            // collect mapping of branch conditions
            var code = body.getCode() instanceof BlockStatement block ? block : null;
            if( code == null )
                return body;

            var allLabels = new LinkedHashSet<String>();
            var allBlocks = new LinkedHashMap<String,BranchCondition>();
            var statements = new ArrayList<Statement>();
            String currentLabel = null;

            for( var stmt : code.getStatements() ) {
                var label = stmt.getStatementLabel();
                if( label != null ) {
                    if( !allLabels.add(label) ) {
                        syntaxError("Branch label already declared: " + label, stmt);
                        break;
                    }
                    currentLabel = label;

                    if( stmt instanceof ExpressionStatement es ) {
                        var block = new BranchCondition(label, es.getExpression());
                        allBlocks.put(label, block);
                    }
                    else {
                        syntaxError("Unexpected statement in label " + label, stmt);
                        break;
                    }
                }
                else if( currentLabel != null ) {
                    var block = allBlocks.get(currentLabel);
                    block.code.add(stmt);
                }
                else {
                    statements.add(stmt);
                }
            }

            if( allBlocks.isEmpty() ) {
                syntaxError("Branch criteria should declare at least one branch", body);
                return body;
            }

            // construct if statement for each branch condition
            IfStatement ifStatement = null;
            IfStatement prevIf = null;
            for( var branch : allBlocks.values() ) {
                var nextIf = ifS(boolX(branch.condition), branchBlock(branch));
                nextIf.addStatementLabel(branch.label);
                if( ifStatement == null )
                    ifStatement = nextIf;
                if( prevIf != null )
                    prevIf.setElseBlock(nextIf);
                prevIf = nextIf;
            }

            statements.add(ifStatement);

            // construct branch criteria
            var main = closureX(body.getParameters(), block(body.getVariableScope(), statements));
            var newTokenBranchDef = createX("nextflow.script.TokenBranchDef", main, list2args(new ArrayList(allLabels)));
            return closureX(null, block(body.getVariableScope(), stmt(newTokenBranchDef)));
        }

        private Statement branchBlock(BranchCondition branch) {
            var choice = branch.label;
            var statements = branch.code;

            // return the closure param by default
            if( statements.isEmpty() ) {
                statements.add(branchReturn(paramX(body.getParameters()), choice));
                return block(body.getVariableScope(), statements);
            }

            // otherwise transform any return statements
            var returns = new BranchReturnCollector().collect(statements);
            if( !returns.isEmpty() ) {
                for( var stmt : returns )
                    stmt.setExpression(branchReturnX(stmt.getExpression(), choice));
                return block(body.getVariableScope(), statements);
            }

            // otherwise transform the last expression statement
            var last = statements.size() - 1;
            if( statements.get(last) instanceof ExpressionStatement es ) {
                statements.set(last, branchReturn(es.getExpression(), choice));
                return block(body.getVariableScope(), statements);
            }

            syntaxError("Unexpected statement in branch condition", statements.get(last));
            return EmptyStatement.INSTANCE;
        }

        private Expression paramX(Parameter[] params) {
            if( params.length == 0 )
                return varX("it");

            if( params.length == 1 )
                return varX(params[0].getName());

            return listX(
                Arrays.stream(params).map(p -> (Expression) varX(p.getName())).toList()
            );
        }

        private Statement branchReturn(Expression value, String choice) {
            return returnS(branchReturnX(value, choice));
        }

        private Expression branchReturnX(Expression value, String choice) {
            return createX("nextflow.script.TokenBranchChoice", value, constX(choice));
        }

        private static class BranchReturnCollector extends CodeVisitorSupport {

            private List<ReturnStatement> returns = new ArrayList<>();

            public List<ReturnStatement> collect(List<Statement> statements) {
                for( var stmt : statements )
                    stmt.visit(this);
                return returns;
            }

            @Override
            public void visitReturnStatement(ReturnStatement node) {
                returns.add(node);
            }
        }

        private static class BranchCondition {
            String label;
            Expression condition;
            List<Statement> code = new ArrayList<>();

            public BranchCondition(String label, Expression condition) {
                this.label = label;
                this.condition = condition;
            }
        }
    }

    class MultiMapTransformer {

        private ClosureExpression body;

        public MultiMapTransformer(ClosureExpression body) {
            this.body = body;
        }

        Expression apply() {
            // assign output variables for each multiMap block
            var code = body.getCode() instanceof BlockStatement block ? block : null;
            if( code == null )
                return body;

            var allLabels = new LinkedHashSet<String>();
            var vars = new LinkedHashSet<String>();
            var statements = new ArrayList<Statement>(code.getStatements().size() * 2);
            List<String> labels = Collections.emptyList();

            for( int i=0; i<code.getStatements().size(); i++ ) {
                var stmt = code.getStatements().get(i);
                var currentLabels = Optional.ofNullable(stmt.getStatementLabels()).orElse(labels);
                if( currentLabels != null )
                    allLabels.addAll(DefaultGroovyMethods.asReversed(currentLabels));
                if( DefaultGroovyMethods.equals(currentLabels, labels) ) {
                    statements.add(stmt);
                    continue;
                }
                if( DefaultGroovyMethods.asBoolean(labels) ) {
                    assignOutputVars(stmt, code.getStatements().get(i - 1), labels, vars, statements);
                    statements.add(stmt);
                }
                else {
                    statements.add(stmt);
                }
                labels = currentLabels;
            }

            if( DefaultGroovyMethods.asBoolean(labels) ) {
                var last = code.getStatements().get(code.getStatements().size() - 1);
                assignOutputVars(last, last, labels, vars, statements);
            }

            if( allLabels.isEmpty() ) {
                syntaxError("multiMap criteria should declare at least two outputs", code);
                return body;
            }

            // construct map entry for each output variable
            statements.add(returnS(mapX(
                vars.stream().map(v -> entryX(constX(v.substring(5)), varX(v))).toList()
            )));

            // construct multiMap criteria
            var main = closureX(body.getParameters(), block(body.getVariableScope(), statements));
            var newTokenMultiMapDef = createX("nextflow.script.TokenMultiMapDef", main, list2args(new ArrayList(allLabels)));
            return closureX(null, block(body.getVariableScope(), stmt(newTokenMultiMapDef)));
        }

        /**
         * Assign the last expression statement in a multiMap block to a variable
         * for each label in the block.
         *
         * This is done to reuse the same result expression for multiple labels.
         *
         * @param current
         * @param previous
         * @param labels
         * @param vars
         * @param statements
         */
        void assignOutputVars(Statement current, Statement previous, List<String> labels, Set<String> vars, List<Statement> statements) {
            VariableExpression target = null;
            for( var label : labels ) {
                var varName = "$out_" + label;
                if( !vars.add(varName) )
                    syntaxError("multiMap label already declared: " + label, current);

                if( !(previous instanceof ExpressionStatement) )
                    syntaxError("multiMap block must end with an expression statement", previous);

                if( target == null ) {
                    var source = ((ExpressionStatement) previous).getExpression();
                    target = varX(varName);
                    statements.add(declS(target, source));
                }
                else {
                    statements.add(declS(varX(varName), target));
                }
            }
        }
    }
}
