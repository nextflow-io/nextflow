/*
 * Copyright 2025, Seqera Labs
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
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import nextflow.script.ast.ASTNodeMarker;
import nextflow.script.ast.AssignmentExpression;
import nextflow.script.ast.ProcessNodeV2;
import nextflow.script.ast.TupleParameter;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.CodeVisitorSupport;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.VariableScope;
import org.codehaus.groovy.ast.expr.BinaryExpression;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.EmptyExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.MapExpression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.BlockStatement;
import org.codehaus.groovy.ast.stmt.EmptyStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.SourceUnit;

import static nextflow.script.ast.ASTUtils.*;
import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Transform a typed process AST node into Groovy AST.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ProcessToGroovyVisitorV2 {

    private SourceUnit sourceUnit;

    private ScriptToGroovyHelper sgh;

    public ProcessToGroovyVisitorV2(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
        this.sgh = new ScriptToGroovyHelper(sourceUnit);
    }

    public Statement transform(ProcessNodeV2 node) {
        visitProcessDirectives(node.directives);
        visitProcessStagers(node.stagers);

        var stagers = node.stagers instanceof BlockStatement block ? block : new BlockStatement();
        visitProcessInputs(node.inputs, stagers);

        var unstagers = new BlockStatement();
        var unstageVisitor = new ProcessUnstageVisitor(unstagers);
        visitProcessUnstagers(node.outputs, unstageVisitor);
        visitProcessUnstagers(node.topics, unstageVisitor);

        if( "script".equals(node.type) )
            node.exec.visit(new TaskCmdXformVisitor(sourceUnit));
        node.stub.visit(new TaskCmdXformVisitor(sourceUnit));

        var body = closureX(block(new VariableScope(), List.of(
            node.directives,
            stagers,
            unstagers,
            processInputs(node.inputs),
            processOutputs(node.outputs),
            processTopics(node.topics),
            processWhen(node.when),
            processStub(node.stub),
            stmt(createX(
                "nextflow.script.BodyDef",
                args(
                    closureX(node.exec),
                    constX(sgh.getSourceText(node.exec)),
                    constX(node.type),
                    sgh.getVariableRefs(node.exec)
                )
            ))
        )));
        return stmt(callThisX("processV2", args(constX(node.getName()), body)));
    }

    private void visitProcessDirectives(Statement directives) {
        asDirectives(directives).forEach((call) -> {
            var arguments = asMethodCallArguments(call);
            if( arguments.size() != 1 )
                return;
            var firstArg = arguments.get(0);
            if( firstArg instanceof ClosureExpression )
                return;
            arguments.set(0, sgh.transformToLazy(firstArg));
        });
    }

    private void visitProcessStagers(Statement directives) {
        asDirectives(directives).forEach((call) -> {
            var arguments = asMethodCallArguments(call).stream()
                .map(arg -> sgh.transformToLazy(arg))
                .toList();
            call.setArguments(args(arguments));
        });
    }

    private void visitProcessInputs(Parameter[] inputs, BlockStatement stagers) {
        for( var param : asFlatParams(inputs) ) {
            if( isPathType(param.getType()) ) {
                var ve = varX(param.getName());
                var stager = stmt(callThisX("stageAs", args(closureX(stmt(ve)))));
                stagers.addStatement(stager);
            }
        }
    }

    private static boolean isPathType(ClassNode cn) {
        if( !cn.isResolved() )
            return false;
        var tn = new TypeNode(cn);
        var type = tn.type;
        if( Path.class.isAssignableFrom(type) ) {
            return true;
        }
        if( Collection.class.isAssignableFrom(type) && tn.genericTypes != null ) {
            var genericType = tn.genericTypes.get(0);
            return Path.class.isAssignableFrom(genericType);
        }
        return false;
    }

    private static class TypeNode {
        final Class type;
        final List<Class> genericTypes;

        public TypeNode(ClassNode cn) {
            this.type = cn.getTypeClass();
            if( cn.isUsingGenerics() ) {
                this.genericTypes = Arrays.stream(cn.getGenericsTypes())
                    .map(el -> el.getType().getTypeClass())
                    .toList();
            }
            else {
                this.genericTypes = null;
            }
        }
    }

    private void visitProcessUnstagers(Statement outputs, ProcessUnstageVisitor visitor) {
        for( var output : asBlockStatements(outputs) )
            visitor.visit(output);
    }

    private static class ProcessUnstageVisitor extends CodeVisitorSupport {

        private int evalCount = 0;

        private int pathCount = 0;

        private BlockStatement unstagers;

        public ProcessUnstageVisitor(BlockStatement unstagers) {
            this.unstagers = unstagers;
        }

        @Override
        public void visitMethodCallExpression(MethodCallExpression node) {
            extractUnstageDirective(node);
            super.visitMethodCallExpression(node);
        }

        private void extractUnstageDirective(MethodCallExpression node) {
            if( !node.isImplicitThis() )
                return;

            var name = node.getMethodAsString();
            var arguments = asMethodCallArguments(node);

            // env(<name>)
            // emit: _unstage_env(<name>)
            if( "env".equals(name) && arguments.size() == 1 ) {
                var key = arguments.get(0);
                var unstager = stmt(callThisX("_unstage_env", args(key)));
                unstagers.addStatement(unstager);

                // rename to _env() to prevent dispatch to equivalent ScriptDsl function
                node.setMethod(constX("_" + name));
            }

            // eval(<cmd>) -> eval(<key>)
            // emit: _unstage_eval(<key>, { <cmd> })
            if( "eval".equals(name) && arguments.size() == 1 ) {
                var key = constX("nxf_out_eval_" + (evalCount++));
                var cmd = arguments.get(0);
                var unstager = stmt(callThisX("_unstage_eval", args(key, closureX(stmt(cmd)))));
                unstagers.addStatement(unstager);
                node.setArguments(args(key));
            }

            // file(<opts>, <pattern>) -> file(<opts>, <key>)
            // files(<opts>, <pattern>) -> files(<opts>, <key>)
            // emit: _unstage_files(<key>, { <pattern> })
            if( "file".equals(name) || "files".equals(name) ) {
                Expression opts;
                Expression pattern;
                if( arguments.size() == 1 ) {
                    opts = new MapExpression();
                    pattern = arguments.get(0);
                }
                else if( arguments.size() == 2 ) {
                    opts = arguments.get(0);
                    pattern = arguments.get(1);
                }
                else {
                    return;
                }

                var key = constX("$path" + (pathCount++));
                var unstager = stmt(callThisX("_unstage_files", args(key, closureX(stmt(pattern)))));
                unstagers.addStatement(unstager);

                // rename to _file() or _files() to prevent dispatch to equivalent ScriptDsl functions
                node.setMethod(constX("_" + name));
                node.setArguments(args(opts, key));
            }
        }
    }

    private Statement processInputs(Parameter[] inputs) {
        var statements = Arrays.stream(inputs)
            .map((input) -> {
                if( input instanceof TupleParameter tp ) {
                    var components = Arrays.stream(tp.components)
                        .map(p -> processInputCtor(p))
                        .toList();
                    var type = input.getType();
                    return stmt(callThisX("_input_", args(listX(components), classX(type))));
                }
                else {
                    var name = input.getName();
                    var type = input.getType();
                    var optional = type.getNodeMetaData(ASTNodeMarker.NULLABLE) != null;
                    return stmt(callThisX("_input_", args(constX(name), classX(type), constX(optional))));
                }
            })
            .toList();
        return block(null, statements);
    }

    private Expression processInputCtor(Parameter param) {
        return createX(
            "nextflow.script.params.v2.ProcessInput",
            args(
                constX(param.getName()),
                classX(param.getType()),
                constX(param.getType().getNodeMetaData(ASTNodeMarker.NULLABLE) != null)
            )
        );
    }

    private Statement processOutputs(Statement outputs) {
        var statements = asBlockStatements(outputs).stream()
            .map(stmt -> ((ExpressionStatement) stmt).getExpression())
            .map((output) -> {
                if( output instanceof VariableExpression target ) {
                    return stmt(callThisX("_output_", args(constX(target.getName()), classX(target.getType()), closureX(stmt(target)))));
                }
                else if( output instanceof AssignmentExpression ae ) {
                    var target = (VariableExpression)ae.getLeftExpression();
                    return stmt(callThisX("_output_", args(constX(target.getName()), classX(target.getType()), closureX(stmt(ae.getRightExpression())))));
                }
                else {
                    return stmt(callThisX("_output_", args(constX("$out"), classX(ClassHelper.dynamicType()), closureX(stmt(output)))));
                }
            })
            .toList();
        return block(null, statements);
    }

    private Statement processTopics(Statement topics) {
        var statements = asBlockStatements(topics).stream()
            .map((stmt) -> {
                var es = (ExpressionStatement) stmt;
                var be = (BinaryExpression) es.getExpression();
                return stmt(callThisX("_topic_", args(closureX(stmt(be.getLeftExpression())), be.getRightExpression())));
            })
            .toList();
        return block(null, statements);
    }

    private Statement processWhen(Expression when) {
        if( when instanceof EmptyExpression )
            return EmptyStatement.INSTANCE;
        return stmt(callThisX("when", createX(
            "nextflow.script.TaskClosure",
            args(
                closureX(null, block(stmt(when))),
                constX(sgh.getSourceText(when))
            )
        )));
    }

    private Statement processStub(Statement stub) {
        if( stub instanceof EmptyStatement )
            return EmptyStatement.INSTANCE;
        return stmt(callThisX("stub", createX(
            "nextflow.script.TaskClosure",
            args(
                closureX(null, stub),
                constX(sgh.getSourceText(stub))
            )
        )));
    }

}
