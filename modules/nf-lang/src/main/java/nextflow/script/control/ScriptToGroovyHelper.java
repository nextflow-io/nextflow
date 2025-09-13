/*
 * Copyright 2025, Seqera Labs
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

import java.util.HashSet;
import java.util.Set;

import org.codehaus.groovy.ast.CodeVisitorSupport;
import org.codehaus.groovy.ast.Variable;
import org.codehaus.groovy.ast.expr.ClosureExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.VariableExpression;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.SourceUnit;

import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Utility functions for ScriptToGroovyVisitor.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ScriptToGroovyHelper {

    private SourceUnit sourceUnit;

    public ScriptToGroovyHelper(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
    }

    /**
     * Transform an expression into a lazy expression by
     * wrapping it in a closure if it references variables.
     *
     * @param node
     */
    public Expression transformToLazy(Expression node)  {
        if( node instanceof ClosureExpression )
            return node;
        var vars = new VariableCollector().collect(node);
        if( !vars.isEmpty() )
            return closureX(stmt(node));
        return node;
    }

    private class VariableCollector extends CodeVisitorSupport {

        private Set<Variable> vars;

        private Set<Variable> declaredParams;

        public Set<Variable> collect(Expression node) {
            vars = new HashSet<>();
            declaredParams = new HashSet<>();
            visit(node);
            return vars;
        }

        @Override
        public void visitClosureExpression(ClosureExpression node) {
            if( node.getParameters() != null ) {
                for( var param : node.getParameters() )
                    declaredParams.add(param);
            }
        }

        @Override
        public void visitVariableExpression(VariableExpression node) {
            var variable = node.getAccessedVariable();
            if( variable != null && !declaredParams.contains(variable) )
                vars.add(variable);
        }
    }

    /**
     * Get the source text for a statement.
     *
     * @param node
     */
    public String getSourceText(Statement node) {
        var builder = new StringBuilder();
        var colx = node.getColumnNumber();
        var colz = node.getLastColumnNumber();
        var first = node.getLineNumber();
        var last = node.getLastLineNumber();
        for( int i = first; i <= last; i++ ) {
            var line = sourceUnit.getSource().getLine(i, null);

            // prepend first-line indent
            if( i == first ) {
                int k = 0;
                while( k < line.length() && line.charAt(k) == ' ' )
                    k++;
                builder.append( line.substring(0, k) );
            }

            var begin = (i == first) ? colx - 1 : 0;
            var end = (i == last) ? colz - 1 : line.length();
            builder.append( line.substring(begin, end) ).append('\n');
        }
        return builder.toString();
    }

    /**
     * Get the source text for an expression.
     *
     * @param node
     */
    public String getSourceText(Expression node) {
        var stm = stmt(node);
        stm.setSourcePosition(node);
        return getSourceText(stm);
    }

}
