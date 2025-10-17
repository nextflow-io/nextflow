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
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;

import static nextflow.script.ast.ASTUtils.*;

/**
 * A typed process definition.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ProcessNodeV2 extends ProcessNode {
    public final Statement directives;
    public final Parameter[] inputs;
    public final Statement stagers;
    public final Statement outputs;
    public final Statement topics;
    public final Expression when;
    public final String type;
    public final Statement exec;
    public final Statement stub;

    public ProcessNodeV2(String name, Statement directives, Parameter[] inputs, Statement stagers, Statement outputs, Statement topics, Expression when, String type, Statement exec, Statement stub) {
        super(name, inputs, dummyReturnType(outputs));
        this.directives = directives;
        this.inputs = inputs;
        this.stagers = stagers;
        this.outputs = outputs;
        this.topics = topics;
        this.when = when;
        this.type = type;
        this.exec = exec;
        this.stub = stub;
    }

    /**
     * Process outputs are represented as a single record, or
     * a value if there is a single output expression.
     *
     * @param block
     */
    private static ClassNode dummyReturnType(Statement block) {
        var outputs = asBlockStatements(block);
        if( outputs.size() == 1 ) {
            var first = outputs.get(0);
            var output = ((ExpressionStatement) first).getExpression();
            if( outputTarget(output) == null )
                return output.getType();
        }
        var cn = new ClassNode(Record.class);
        outputs.stream()
            .map(stmt -> ((ExpressionStatement) stmt).getExpression())
            .map(output -> outputTarget(output))
            .filter(target -> target != null)
            .forEach((target) -> {
                var fn = new FieldNode(target.getName(), Modifier.PUBLIC, target.getType(), cn, null);
                fn.setDeclaringClass(cn);
                cn.addField(fn);
            });
        return cn;
    }

    private static VariableExpression outputTarget(Expression output) {
        if( output instanceof VariableExpression ve ) {
            return ve;
        }
        if( output instanceof AssignmentExpression ae ) {
            return (VariableExpression)ae.getLeftExpression();
        }
        return null;
    }
}
