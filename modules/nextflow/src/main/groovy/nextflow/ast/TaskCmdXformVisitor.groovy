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

import static org.codehaus.groovy.ast.tools.GeneralUtils.*

import groovy.transform.CompileStatic
import org.codehaus.groovy.ast.ClassCodeVisitorSupport
import org.codehaus.groovy.ast.expr.Expression
import org.codehaus.groovy.ast.expr.GStringExpression
import org.codehaus.groovy.ast.tools.GeneralUtils
import org.codehaus.groovy.control.SourceUnit
/**
 * Implements an AST xform visitor to escape and manipulate task command
 * special value
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TaskCmdXformVisitor extends ClassCodeVisitorSupport {

    private SourceUnit unit

    TaskCmdXformVisitor(SourceUnit unit) { this.unit = unit }

    @Override
    protected SourceUnit getSourceUnit() { unit }

    /**
     * Intercepts Gstring commands to apply special formatting
     * such escaping blanks in path
     */
    @Override
    void visitGStringExpression(GStringExpression expr) {
        for(int i=0; i<expr.values.size(); i++)
            expr.values[i] = interceptX(expr.values[i])
        // proceed with the gstring invocation
        super.visitGStringExpression(expr)
    }

    /**
     * {@link LangHelpers#applyPathEscapeAware(java.lang.Object)}
     */
    static protected Expression interceptX(Expression expr) {
        GeneralUtils.callX(
                classX(LangHelpers), 'applyPathEscapeAware',args(expr))
    }
}
