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
package nextflow.script.ast;

import java.lang.reflect.Modifier;

import nextflow.script.types.Record;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.FieldNode;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.EmptyStatement;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;

import static nextflow.script.ast.ASTUtils.*;

/**
 * A workflow definition.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class WorkflowNode extends MethodNode {
    public final Statement main;
    public final Statement emits;
    public final Statement publishers;
    public final Statement onComplete;
    public final Statement onError;

    public WorkflowNode(String name, Parameter[] takes, Statement main, Statement emits, Statement publishers, Statement onComplete, Statement onError) {
        super(name, 0, dummyReturnType(emits), takes, ClassNode.EMPTY_ARRAY, EmptyStatement.INSTANCE);
        this.main = main;
        this.emits = emits;
        this.publishers = publishers;
        this.onComplete = onComplete;
        this.onError = onError;
    }

    public WorkflowNode(String name, Statement main) {
        this(name, Parameter.EMPTY_ARRAY, main, EmptyStatement.INSTANCE, EmptyStatement.INSTANCE, EmptyStatement.INSTANCE, EmptyStatement.INSTANCE);
    }

    public boolean isEntry() {
        return getName() == null;
    }

    public boolean isCodeSnippet() {
        return getLineNumber() == -1;
    }

    private static ClassNode dummyReturnType(Statement block) {
        var emits = asBlockStatements(block);
        if( emits.size() == 1 ) {
            var first = emits.get(0);
            var emit = ((ExpressionStatement) first).getExpression();
            if( emitTarget(emit) == null )
                return emit.getType();
        }
        var cn = new ClassNode(Record.class);
        emits.stream()
            .map(stmt -> ((ExpressionStatement) stmt).getExpression())
            .map(emit -> emitTarget(emit))
            .filter(target -> target != null)
            .forEach((target) -> {
                var fn = new FieldNode(target.getName(), Modifier.PUBLIC, target.getType(), cn, null);
                fn.setDeclaringClass(cn);
                cn.addField(fn);
            });
        return cn;
    }

    private static VariableExpression emitTarget(Expression emit) {
        if( emit instanceof VariableExpression ve ) {
            return ve;
        }
        if( emit instanceof AssignmentExpression ae ) {
            return (VariableExpression)ae.getLeftExpression();
        }
        return null;
    }
}
