/*
 * Copyright 2024, Seqera Labs
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
package nextflow.config.parser.v2;

import java.util.ArrayList;
import java.util.stream.Collectors;

import nextflow.config.ast.ConfigAppendNode;
import nextflow.config.ast.ConfigAssignNode;
import nextflow.config.ast.ConfigBlockNode;
import nextflow.config.ast.ConfigIncludeNode;
import nextflow.config.ast.ConfigNode;
import nextflow.config.ast.ConfigVisitorSupport;
import org.codehaus.groovy.ast.VariableScope;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.stmt.ReturnStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.SourceUnit;

import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Visitor to convert a Nextflow config AST into a
 * Groovy AST which is executed against {@link ConfigDsl}.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ConfigToGroovyVisitor extends ConfigVisitorSupport {

    private SourceUnit sourceUnit;

    private ConfigNode moduleNode;

    public ConfigToGroovyVisitor(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
        this.moduleNode = (ConfigNode) sourceUnit.getAST();
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    public void visit() {
        if( moduleNode == null )
            return;
        super.visit(moduleNode);
        if( moduleNode.isEmpty() )
            moduleNode.addStatement(ReturnStatement.RETURN_NULL_OR_VOID);
    }

    @Override
    public void visitConfigAssign(ConfigAssignNode node) {
        moduleNode.addStatement(transformConfigAssign(node));
    }

    protected Statement transformConfigAssign(ConfigAssignNode node) {
        if( node instanceof ConfigAppendNode ) {
            var method = node.names.get(0);
            return stmt(callThisX(method, args(node.value)));
        }
        var names = listX(
            node.names.stream()
                .map(name -> (Expression) constX(name))
                .collect(Collectors.toList())
        );
        return stmt(callThisX("assign", args(names, node.value)));
    }

    @Override
    public void visitConfigBlock(ConfigBlockNode node) {
        moduleNode.addStatement(transformConfigBlock(node));
    }

    protected Statement transformConfigBlock(ConfigBlockNode node) {
        var statements = new ArrayList<Statement>();
        for( var stmt : node.statements ) {
            if( stmt instanceof ConfigAssignNode can )
                statements.add(transformConfigAssign(can));
            else if( stmt instanceof ConfigBlockNode cbn )
                statements.add(transformConfigBlock(cbn));
            else if( stmt instanceof ConfigIncludeNode cin )
                statements.add(transformConfigInclude(cin));
        }
        var code = block(new VariableScope(), statements);
        var kind = node.kind != null ? node.kind : "block";
        return stmt(callThisX(kind, args(constX(node.name), closureX(code))));
    }

    @Override
    public void visitConfigInclude(ConfigIncludeNode node) {
        moduleNode.addStatement(transformConfigInclude(node));
    }

    protected Statement transformConfigInclude(ConfigIncludeNode node) {
        return stmt(callThisX("includeConfig", args(node.source)));
    }

}
