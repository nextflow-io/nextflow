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

package nextflow.config.parser.v2;

import nextflow.config.ConfigClosurePlaceholder;
import org.codehaus.groovy.ast.ClassCodeVisitorSupport;
import org.codehaus.groovy.ast.ClassNode;
import org.codehaus.groovy.ast.expr.ArgumentListExpression;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.MapExpression;
import org.codehaus.groovy.ast.expr.MethodCallExpression;
import org.codehaus.groovy.control.SourceUnit;

import static org.codehaus.groovy.ast.tools.GeneralUtils.*;
/**
 * AST transformation to render closure source text
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ClosureToStringVisitor extends ClassCodeVisitorSupport {

    protected SourceUnit sourceUnit;

    public ClosureToStringVisitor(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    @Override
    public void visitMethodCallExpression(MethodCallExpression node) {
        var name = node.getMethodAsString();
        if( !"assign".equals(name) ) {
            super.visitMethodCallExpression(node);
            return;
        }

        var arguments = (ArgumentListExpression) node.getArguments();
        if( arguments.getExpressions().size() != 2 )
            return;

        var secondArg = arguments.getExpression(1);
        if( secondArg instanceof MapExpression me ) {
            for( var entry : me.getMapEntryExpressions() ) {
                if( entry.getValueExpression() instanceof ClosureExpression ce )
                    entry.setValueExpression(wrapExprAsPlaceholder(ce));
            }
        }
        else if( secondArg instanceof ClosureExpression ce ) {
            node.setArguments(args(arguments.getExpression(0), wrapExprAsPlaceholder(ce)));
        }
    }

    protected Expression wrapExprAsPlaceholder(Expression node) {
        var type = new ClassNode(ConfigClosurePlaceholder.class);
        var text = constX(getSourceText(node));
        return ctorX(type, args(text));
    }

    protected String getSourceText(Expression node) {
        var builder = new StringBuilder();
        var colBegin = Math.max(node.getColumnNumber()-1, 0);
        var colEnd = Math.max(node.getLastColumnNumber()-1, 0);
        var lineFirst = node.getLineNumber();
        var lineLast = node.getLastLineNumber();

        for( int i = lineFirst; i <= lineLast; i++ ) {
            var line = sourceUnit.getSource().getLine(i, null);
            if( i == lineFirst ) {
                var str = i == lineLast ? line.substring(colBegin,colEnd) : line.substring(colBegin);
                builder.append(str);
            }
            else {
                var str = i == lineLast ? line.substring(0, colEnd) : line;
                builder.append('\n');
                builder.append(str);
            }
        }
        return builder.toString();
    }

}
