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

import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.MethodNode;

/**
 * An included process, workflow, or function.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class IncludeModuleNode extends ASTNode {
    public final String name;
    public final String alias;

    public IncludeModuleNode(String name, String alias) {
        this.name = name;
        this.alias = alias;
    }

    public IncludeModuleNode(String name) {
        this(name, null);
    }

    public String getNameOrAlias() {
        return alias != null ? alias : name;
    }

    private MethodNode target;

    public void setTarget(MethodNode target) {
        this.target = target;
    }

    public MethodNode getTarget() {
        return target;
    }
}
