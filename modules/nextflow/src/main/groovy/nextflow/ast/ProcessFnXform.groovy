/*
 * Copyright 2013-2023, Seqera Labs
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

import static org.codehaus.groovy.ast.tools.GeneralUtils.*

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.processor.TaskConfig
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.ClassCodeVisitorSupport
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.MethodNode
import org.codehaus.groovy.ast.Parameter
import org.codehaus.groovy.ast.expr.ArgumentListExpression
import org.codehaus.groovy.ast.expr.BinaryExpression
import org.codehaus.groovy.ast.expr.ClosureExpression
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.ListExpression
import org.codehaus.groovy.ast.expr.MethodCallExpression
import org.codehaus.groovy.ast.expr.UnaryMinusExpression
import org.codehaus.groovy.ast.expr.VariableExpression
import org.codehaus.groovy.ast.stmt.BlockStatement
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.ast.stmt.Statement
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.syntax.SyntaxException
import org.codehaus.groovy.syntax.Types
/**
 * Implements syntax transformations for process functions.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ProcessFnXform extends ClassCodeVisitorSupport {

    private final SourceUnit unit

    ProcessFnXform(SourceUnit unit) {
        this.unit = unit
    }

    @Override
    protected SourceUnit getSourceUnit() { unit }

    @Override
    void visitMethod(MethodNode method) {
        final annotation = method.getAnnotations()
                .find(a -> a.getClassNode().getName() == 'ProcessFn')

        if( annotation )
            transform(method, annotation.getMembers())
    }

    /**
     * Transform a ProcessFn annotation to be semantically valid.
     *
     * @param method
     * @param opts
     */
    protected void transform(MethodNode method, Map<String,Expression> opts) {
        // fix directives
        final directives = opts['directives']
        if( directives != null && directives instanceof ClosureExpression ) {
            final block = (BlockStatement)directives.getCode()
            for( Statement stmt : block.getStatements() )
                fixDirectiveWithNegativeValue(stmt)
        }

        // fix outputs
        final outputs = opts['outputs']
        if( outputs != null && outputs instanceof ClosureExpression ) {
            final block = (BlockStatement)outputs.getCode()
            for( Statement stmt : block.getStatements() )
                fixOutputMethod((ExpressionStatement)stmt)
        }

        // insert `task` method parameter
        final params = method.getParameters() as List<Parameter>
        params.push(new Parameter(new ClassNode(TaskConfig), 'task'))
        method.setParameters(params as Parameter[])

        // TODO: append stub source

        // append method params
        opts.put( 'params', closureX( block( new ExpressionStatement(
            new ListExpression(
                params.collect(p -> (Expression)constX(p.getName()))
            )
        ) ) ) )

        // append script source
        opts.put( 'source', constX( getSource(method.getCode()) ) )
    }

    /**
     * Fix directives with a single argument with a negative value,
     * since it will be parsed as a subtract expression if there are
     * no parentheses.
     *
     * @param stmt
     */
    protected void fixDirectiveWithNegativeValue(Statement stmt) {
        // -- check for binary subtract expression
        if( stmt !instanceof ExpressionStatement )
            return
        def expr = ((ExpressionStatement)stmt).getExpression()
        if( expr !instanceof BinaryExpression )
            return
        def binary = (BinaryExpression)expr
        if( binary.leftExpression !instanceof VariableExpression )
            return
        if( binary.operation.type != Types.MINUS )
            return

        // -- transform binary expression `NAME - ARG` into method call `NAME(-ARG)`
        def name = ((VariableExpression)binary.leftExpression).name
        def arg = (Expression)new UnaryMinusExpression(binary.rightExpression)

        ((ExpressionStatement)stmt).setExpression( new MethodCallExpression(
            VariableExpression.THIS_EXPRESSION,
            name,
            new ArgumentListExpression(arg)
        ) )
    }

    private static final VALID_OUTPUT_METHODS = ['val','env','file','path','stdout','tuple']

    /**
     * Fix output method calls.
     *
     * @param stmt
     */
    protected void fixOutputMethod(ExpressionStatement stmt) {
        final methodCall = (MethodCallExpression)stmt.getExpression()
        final name = methodCall.getMethodAsString()
        final args = (ArgumentListExpression)methodCall.getArguments()

        if( name !in VALID_OUTPUT_METHODS )
            syntaxError(stmt, "Invalid output method '${name}'")

        methodCall.setMethod( constX('_out_' + name) )
    }

    private String getSource(ASTNode node) {
        final buffer = new StringBuilder()
        final colx = node.getColumnNumber()
        final colz = node.getLastColumnNumber()
        final first = node.getLineNumber()
        final last = node.getLastLineNumber()
        for( int i=first; i<=last; i++ ) {
            def line = unit.source.getLine(i, null)
            if( i==last )
                line = line.substring(0,colz-1)
            if( i==first )
                line = line.substring(colx-1)
            buffer.append(line) .append('\n')
        }

        return buffer.toString()
    }

    protected void syntaxError(ASTNode node, String message) {
        unit.addError( new SyntaxException(message, node.lineNumber, node.columnNumber) )
    }

}
