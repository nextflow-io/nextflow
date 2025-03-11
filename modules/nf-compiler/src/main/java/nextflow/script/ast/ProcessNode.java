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

import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.stmt.EmptyStatement;
import org.codehaus.groovy.ast.stmt.Statement;

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ProcessNode extends MethodNode {
    public final Statement directives;
    public final Statement inputs;
    public final Statement outputs;
    public final Expression when;
    public final String type;
    public final Statement exec;
    public final Statement stub;

    public ProcessNode(String name, Statement directives, Statement inputs, Statement outputs, Expression when, String type, Statement exec, Statement stub) {
        super(name, 0, null, Parameter.EMPTY_ARRAY, ClassNode.EMPTY_ARRAY, EmptyStatement.INSTANCE);
        this.directives = directives;
        this.inputs = inputs;
        this.outputs = outputs;
        this.when = when;
        this.type = type;
        this.exec = exec;
        this.stub = stub;
    }
}
