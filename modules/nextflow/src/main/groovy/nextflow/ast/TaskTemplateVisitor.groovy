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
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.ConstructorNode
import org.codehaus.groovy.ast.expr.ArgumentListExpression
import org.codehaus.groovy.ast.expr.ConstantExpression
import org.codehaus.groovy.ast.expr.GStringExpression
import org.codehaus.groovy.ast.expr.ListExpression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.ast.stmt.Statement
import org.codehaus.groovy.control.SourceUnit
/**
 * Inject string var names in the script binding object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TaskTemplateVisitor extends VariableVisitor {

    int withinGString

    TaskTemplateVisitor(SourceUnit unit) {
        super(unit)
    }

    @Override
    protected boolean scopeEnabled() {
        return withinGString>0
    }

    @Override
    void visitGStringExpression(GStringExpression expression) {
        try {
            withinGString++
            super.visitGStringExpression(expression)
        }
        finally {
            withinGString--
        }
    }

    @Override
    protected void visitObjectInitializerStatements(ClassNode node) {
        injectMetadata(node)
        super.visitObjectInitializerStatements(node)
    }

    protected void injectMetadata(ClassNode node) {
        // inject the invocation of a custom method
        // as last statement of the constructor invocation
        for( ConstructorNode constructor : node.getDeclaredConstructors() ) {
            def code = constructor.getCode()
            if( code instanceof BlockStatement ) {
                code.addStatement(makeSetVariableNamesStm())
            }
            else if( code instanceof ExpressionStatement ) {
                def expr = code
                def block = new BlockStatement()
                block.addStatement(expr)
                block.addStatement(makeSetVariableNamesStm())
                constructor.setCode(block)
            }
            else
                throw new IllegalStateException("Invalid constructor expression: $code")
        }
    }

    protected Statement makeSetVariableNamesStm() {
        final names = new ListExpression()
        for( String it: getAllVariableNames() ) {
            names.addExpression(new ConstantExpression(it))
        }


        final varArgs = new ArgumentListExpression()
        varArgs.addExpression( new ConstantExpression('__$$_template_vars') )
        varArgs.addExpression( names )

        // some magic code
        // this generates the invocation of the method:
        //   this.getBinding().setVariable('__$$_template_vars', [variable names])
        final thiz = new VariableExpression('this')
        final bind = new MethodCallExpression( thiz, 'getBinding', ArgumentListExpression.EMPTY_ARGUMENTS )
        final call = new MethodCallExpression( bind, 'setVariable', varArgs )
        final stm = new ExpressionStatement(call)
        return stm
    }
}
