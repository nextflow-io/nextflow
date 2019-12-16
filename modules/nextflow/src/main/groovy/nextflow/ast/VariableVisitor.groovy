/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.ast

import groovy.transform.CompileStatic
import nextflow.script.TokenValRef
import org.codehaus.groovy.ast.ClassCodeVisitorSupport
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.ConstantExpression
import org.codehaus.groovy.ast.expr.DeclarationExpression
import org.codehaus.groovy.ast.expr.PropertyExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.syntax.SyntaxException
/**
 * Visit a closure and collect all referenced variable names
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class VariableVisitor extends ClassCodeVisitorSupport {

    final Map<String, TokenValRef> fAllVariables = [:]

    final Set<String> localDef = []

    final SourceUnit sourceUnit

    private boolean declaration

    private int deep

    VariableVisitor( SourceUnit unit ) {
        this.sourceUnit = unit
    }

    protected boolean scopeEnabled() { true }

    protected boolean isNormalized(PropertyExpression expr) {
        if( !(expr.getProperty() instanceof ConstantExpression) )
            return false

        def target = expr.getObjectExpression()
        while( target instanceof PropertyExpression) {
            target = (target as PropertyExpression).getObjectExpression()
        }

        return target instanceof VariableExpression
    }

    @Override
    void visitClosureExpression(ClosureExpression expression) {
        if( deep++ == 0 )
            super.visitClosureExpression(expression)
    }

    @Override
    void visitDeclarationExpression(DeclarationExpression expr) {
        declaration = true
        try {
            super.visitDeclarationExpression(expr)
        }
        finally {
            declaration = false
        }
    }

    @Override
    void visitPropertyExpression(PropertyExpression expr) {

        if( isNormalized(expr)) {
            final name = expr.text.replace('?','')
            final line = expr.lineNumber
            final coln = expr.columnNumber

            if( !name.startsWith('this.') && !fAllVariables.containsKey(name) && scopeEnabled() ) {
                fAllVariables[name] = new TokenValRef(name,line,coln)
            }
        }
        else
            super.visitPropertyExpression(expr)

    }

    @Override
    void visitVariableExpression(VariableExpression var) {
        final name = var.name
        final line = var.lineNumber
        final coln = var.columnNumber

        if( name == 'this' )
            return

        if( declaration ) {
            if( fAllVariables.containsKey(name) )
                sourceUnit.addError( new SyntaxException("Variable `$name` already defined in the process scope", line, coln))
            else
                localDef.add(name)
        }

        // Note: variable declared in the process scope are not added
        // to the set of referenced variables. Only global ones are tracked
        else if( !localDef.contains(name) && !fAllVariables.containsKey(name) && scopeEnabled() ) {
            fAllVariables[name] = new TokenValRef(name,line,coln)
        }
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit
    }

    /**
     * @return The set of all variables referenced in the script.
     * NOTE: it includes properties in the form {@code object.propertyName}
     */
    Set<TokenValRef> getAllVariables() {
        new HashSet<TokenValRef>(fAllVariables.values())
    }

    Set<String> getAllVariableNames() {
        def all = getAllVariables()
        def result = new HashSet(all.size())
        for( TokenValRef ref : all )
            result.add(ref.name)
        return result
    }
}
