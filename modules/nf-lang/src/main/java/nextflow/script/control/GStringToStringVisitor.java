/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.script.control;

import org.codehaus.groovy.ast.ClassCodeExpressionTransformer;
import org.codehaus.groovy.ast.expr.ArgumentListExpression;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.GStringExpression;
import org.codehaus.groovy.control.SourceUnit;

import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Coerce a GString into a String.
 *
 * from
 *   "${foo} ${bar}"
 * to
 *   "${foo} ${bar}".toString()
 *
 * This enables equality checks between GStrings and Strings,
 * e.g. `"${'hello'}" == 'hello'`.
 *
 * @author Ben Sherman <bentshermman@gmail.com>
 */
public class GStringToStringVisitor extends ClassCodeExpressionTransformer {

    private SourceUnit sourceUnit;

    public GStringToStringVisitor(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    @Override
    public Expression transform(Expression node) {
        if( node instanceof ClosureExpression ce )
            super.visitClosureExpression(ce);

        if( node instanceof GStringExpression gse )
            return transformToString(gse);

        return super.transform(node);
    }

    private Expression transformToString(GStringExpression node) {
        return callX(node, "toString", new ArgumentListExpression());
    }

}
