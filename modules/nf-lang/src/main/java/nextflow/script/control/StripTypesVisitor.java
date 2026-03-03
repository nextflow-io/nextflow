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

package nextflow.script.control;

import java.util.List;

import nextflow.script.types.Record;
import nextflow.script.types.Tuple;
import nextflow.script.types.Value;
import org.codehaus.groovy.ast.ClassCodeExpressionTransformer;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.expr.CastExpression;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.DeclarationExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.control.SourceUnit;

/**
 * Strip type annotations that are used by the Nextflow type checker
 * but not supported by the Groovy runtime.
 *
 * For example, Nextflow allows the Tuple type to be specified
 * with variable type arguments, which is not supported by the JVM,
 * so these type annotations must be removed.
 *
 * @author Ben Sherman <bentshermman@gmail.com>
 */
public class StripTypesVisitor extends ClassCodeExpressionTransformer {

    private static final List<ClassNode> STRIP_TYPES = List.of(
        ClassHelper.makeWithoutCaching("nextflow.Channel"),
        ClassHelper.makeCached(Record.class),
        ClassHelper.makeCached(Tuple.class),
        ClassHelper.makeCached(Value.class)
    );

    private SourceUnit sourceUnit;

    public StripTypesVisitor(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    @Override
    public Expression transform(Expression node) {
        if( node instanceof CastExpression ce ) {
            return stripTypeAnnotation(ce);
        }

        if( node instanceof ClosureExpression ce ) {
            ce.visit(this);
            return ce;
        }

        if( node instanceof DeclarationExpression de ) {
            stripTypeAnnotation(de);
        }

        return super.transform(node);
    }

    private Expression stripTypeAnnotation(CastExpression node) {
        return STRIP_TYPES.contains(node.getType())
            ? node.getExpression()
            : node;
    }

    private Expression stripTypeAnnotation(DeclarationExpression node) {
        if( node.getLeftExpression() instanceof VariableExpression ve ) {
            if( STRIP_TYPES.contains(ve.getType()) )
                node.setLeftExpression(new VariableExpression(ve.getName()));
        }
        return node;
    }

}
