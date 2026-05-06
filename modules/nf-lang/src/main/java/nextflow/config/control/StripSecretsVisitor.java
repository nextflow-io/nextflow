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

import org.codehaus.groovy.ast.ClassCodeExpressionTransformer;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.ConstantExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.ast.expr.PropertyExpression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.control.SourceUnit;

import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Replace secret expressions with a string literal.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class StripSecretsVisitor extends ClassCodeExpressionTransformer {

    private SourceUnit sourceUnit;

    private boolean isIncludeConfigArgument;

    public StripSecretsVisitor(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    @Override
    public Expression transform(Expression node) {
        if( node instanceof ClosureExpression ce ) {
            ce.visit(this);
            return ce;
        }

        if( node instanceof MethodCallExpression mce ) {
            return transformMethodCall(mce);
        }

        if( node instanceof PropertyExpression pe ) {
            return transformProperty(pe);
        }

        return super.transform(node);
    }

    /**
     * Don't obfuscate secret references in a config include, since
     * it would break the config inclusion and isn't needed anyway
     * (config include source is not preserved in config map).
     *
     * @param node
     */
    private Expression transformMethodCall(MethodCallExpression node) {
        if( "includeConfig".equals(node.getMethodAsString()) ) {
            isIncludeConfigArgument = true;
            try {
                return super.transform(node);
            }
            finally {
                isIncludeConfigArgument = false;
            }
        }
        else {
            return super.transform(node);
        }
    }

    /**
     * Replace any reference to a secret with a string literal
     * in order to not dislose the secret value when printin the config.
     *
     * @param node
     */
    private Expression transformProperty(PropertyExpression node) {
        if( isSecretProperty(node) && !isIncludeConfigArgument )
            return constX("secrets." + node.getPropertyAsString());
        else
            return super.transform(node);
    }

    private boolean isSecretProperty(PropertyExpression node) {
        return node.getObjectExpression() instanceof VariableExpression ve
            && "secrets".equals(ve.getText())
            && node.getProperty() instanceof ConstantExpression;
    }

}
