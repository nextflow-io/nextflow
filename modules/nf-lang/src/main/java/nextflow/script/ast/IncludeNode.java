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

import java.util.List;

import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.ast.expr.ConstantExpression;

/**
 * An include declaration.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class IncludeNode extends ASTNode {
    public final ConstantExpression source;
    public final List<IncludeModuleNode> modules;

    public IncludeNode(ConstantExpression source, List<IncludeModuleNode> modules) {
        this.source = source;
        this.modules = modules;
    }
}
