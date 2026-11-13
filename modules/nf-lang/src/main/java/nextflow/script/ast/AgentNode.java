/*
 * Copyright 2013-2026, Seqera Labs
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

import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.stmt.EmptyStatement;
import org.codehaus.groovy.ast.stmt.Statement;

/**
 * AST node for an agent definition.
 *
 * Mirrors the shape of {@link ProcessNodeV2} but with a simpler body:
 * directives, typed inputs, typed outputs, and a prompt block.
 * No script/exec/stub/when/stage/topic — execution is delegated to the
 * agent runner (see nf-agent plugin).
 */
public class AgentNode extends MethodNode {

    public final Statement directives;
    public final Parameter[] inputs;
    public final Statement outputs;
    public final Statement prompt;

    public AgentNode(String name, Statement directives, Parameter[] inputs, Statement outputs, Statement prompt) {
        super(name, 0, ClassHelper.OBJECT_TYPE, inputs, ClassNode.EMPTY_ARRAY, EmptyStatement.INSTANCE);
        this.directives = directives;
        this.inputs = inputs;
        this.outputs = outputs;
        this.prompt = prompt;
    }
}
