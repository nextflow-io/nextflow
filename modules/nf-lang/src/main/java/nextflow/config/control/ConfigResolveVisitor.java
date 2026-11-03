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
package nextflow.config.control;

import java.util.Collections;
import java.util.List;
import java.util.regex.Pattern;

import nextflow.config.ast.ConfigAssignNode;
import nextflow.config.ast.ConfigIncludeNode;
import nextflow.config.ast.ConfigNode;
import nextflow.config.ast.ConfigVisitorSupport;
import nextflow.script.control.ResolveVisitor;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.DynamicVariable;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.control.CompilationUnit;
import org.codehaus.groovy.control.SourceUnit;

/**
 * Resolve variable names, function names, and type names in
 * a config file.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ConfigResolveVisitor extends ConfigVisitorSupport {

    private SourceUnit sourceUnit;

    private ResolveVisitor resolver;

    public ConfigResolveVisitor(SourceUnit sourceUnit, CompilationUnit compilationUnit, List<ClassNode> defaultImports) {
        this.sourceUnit = sourceUnit;
        this.resolver = new ResolveVisitor(sourceUnit, compilationUnit, defaultImports, Collections.emptyList());
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    public void visit() {
        var moduleNode = sourceUnit.getAST();
        if( moduleNode instanceof ConfigNode cn ) {
            // initialize variable scopes
            new VariableScopeVisitor(sourceUnit).visit();

            // resolve type names
            super.visit(cn);

            // report errors for any unresolved variable references
            new DynamicVariablesVisitor().visit(cn);
        }
    }

    @Override
    public void visitConfigAssign(ConfigAssignNode node) {
        node.value = resolver.transform(node.value);
    }

    @Override
    public void visitConfigInclude(ConfigIncludeNode node) {
        node.source = resolver.transform(node.source);
    }

    private class DynamicVariablesVisitor extends ConfigVisitorSupport {

        private static final Pattern ENV_VAR_NAME = Pattern.compile("[A-Z_]+[A-Z0-9_]*");

        @Override
        protected SourceUnit getSourceUnit() {
            return sourceUnit;
        }

        @Override
        public void visitVariableExpression(VariableExpression node) {
            var variable = node.getAccessedVariable();
            if( variable instanceof DynamicVariable ) {
                var message = "`" + node.getName() + "` is not defined";
                if( ENV_VAR_NAME.matcher(variable.getName()).matches() )
                    message += " (hint: use `env('...')` to access environment variable)";
                resolver.addError(message, node);
            }
        }
    }

}
