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
import java.util.Optional;

import nextflow.script.types.Channel;
import nextflow.script.types.NamedTuple;
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
    public final Statement takes;
    public final Statement main;
    public final Statement emits;
    public final Statement publishers;

    public WorkflowNode(String name, Statement takes, Statement main, Statement emits, Statement publishers) {
        super(name, 0, dummyReturnType(emits), Parameter.EMPTY_ARRAY, ClassNode.EMPTY_ARRAY, EmptyStatement.INSTANCE);
        this.takes = takes;
        this.main = main;
        this.emits = emits;
        this.publishers = publishers;
    }

    public boolean isEntry() {
        return getName() == null;
    }

    private static ClassNode dummyReturnType(Statement emits) {
        var cn = new ClassNode(NamedTuple.class);
        asBlockStatements(emits).stream()
            .map(stmt -> ((ExpressionStatement) stmt).getExpression())
            .map(emit -> emitName(emit))
            .filter(name -> name != null)
            .forEach((name) -> {
                var type = ClassHelper.dynamicType();
                var fn = new FieldNode(name, Modifier.PUBLIC, type, cn, null);
                fn.setDeclaringClass(cn);
                cn.addField(fn);
            });
        return cn;
    }

    private static String emitName(Expression emit) {
        if( emit instanceof VariableExpression ve ) {
            return ve.getName();
        }
        else if( emit instanceof AssignmentExpression ae ) {
            var left = (VariableExpression)ae.getLeftExpression();
            return left.getName();
        }
        return null;
    }
}
