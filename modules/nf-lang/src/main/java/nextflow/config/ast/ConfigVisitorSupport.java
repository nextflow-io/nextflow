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
package nextflow.config.ast;

import org.codehaus.groovy.ast.ClassCodeVisitorSupport;
import org.codehaus.groovy.ast.expr.MethodCallExpression;

public abstract class ConfigVisitorSupport extends ClassCodeVisitorSupport implements ConfigVisitor {

    //--------------------------------------------------------------------------
    // config statements

    @Override
    public void visit(ConfigNode node) {
        for( var statement : node.getConfigStatements() ) {
            visit(statement);
        }
    }

    @Override
    public void visit(ConfigStatement node) {
        node.visit(this);
    }

    @Override
    public void visitConfigApplyBlock(ConfigApplyBlockNode node) {
        for( var statement : node.statements ) {
            visitConfigApply(statement);
        }
    }

    public void visitConfigApply(ConfigApplyNode node) {
    }

    @Override
    public void visitConfigAssign(ConfigAssignNode node) {
        visit(node.value);
    }

    @Override
    public void visitConfigBlock(ConfigBlockNode node) {
        for( var statement : node.statements ) {
            visit(statement);
        }
    }

    @Override
    public void visitConfigInclude(ConfigIncludeNode node) {
        visit(node.source);
    }

    @Override
    public void visitConfigIncomplete(ConfigIncompleteNode node) {
    }

    //--------------------------------------------------------------------------
    // expressions

    @Override
    public void visitMethodCallExpression(MethodCallExpression node) {
        if( !node.isImplicitThis() )
            node.getObjectExpression().visit(this);
        node.getMethod().visit(this);
        node.getArguments().visit(this);
    }

}
