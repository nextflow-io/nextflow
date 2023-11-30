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
import org.codehaus.groovy.ast.ASTNode
import org.codehaus.groovy.ast.AnnotationNode
import org.codehaus.groovy.ast.ClassCodeVisitorSupport
import org.codehaus.groovy.ast.MethodNode
import org.codehaus.groovy.ast.Parameter
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.ListExpression
import org.codehaus.groovy.ast.stmt.ExpressionStatement
import org.codehaus.groovy.control.SourceUnit
/**
 * Implements syntax transformations for workflow functions.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class WorkflowFnXform extends ClassCodeVisitorSupport {

    private final SourceUnit unit

    WorkflowFnXform(SourceUnit unit) {
        this.unit = unit
    }

    @Override
    protected SourceUnit getSourceUnit() { unit }

    @Override
    void visitMethod(MethodNode method) {
        final annotation = method.getAnnotations()
                .find(a -> a.getClassNode().getName() == 'WorkflowFn')

        if( annotation )
            transform(method, annotation)
    }

    protected void transform(MethodNode method, AnnotationNode annotation) {
        // append method params
        final params = method.getParameters() as List<Parameter>
        annotation.addMember( 'params', closureX( block( new ExpressionStatement(
            new ListExpression(
                params.collect(p -> (Expression)constX(p.getName()))
            )
        ) ) ) )

        // append workflow source
        annotation.addMember( 'source', constX( getSource(method.getCode()) ) )
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

}
