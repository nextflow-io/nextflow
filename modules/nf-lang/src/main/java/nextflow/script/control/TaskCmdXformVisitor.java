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

import org.codehaus.groovy.ast.ClassCodeVisitorSupport;
import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.expr.GStringExpression;
import org.codehaus.groovy.control.SourceUnit;

import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Transform a process script to use Bash-aware path escaping.
 *
 * This way, the user can reference files in a Bash script
 * without needing to e.g. escape spaces in filenames.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class TaskCmdXformVisitor extends ClassCodeVisitorSupport {

    private SourceUnit sourceUnit;

    public TaskCmdXformVisitor(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
    }

    @Override
    protected SourceUnit getSourceUnit() {
        return sourceUnit;
    }

    @Override
    public void visitGStringExpression(GStringExpression node) {
        var values = node.getValues();
        for( int i = 0; i < values.size(); i++ )
            values.set(i, applyEscape(values.get(i)));
        super.visitGStringExpression(node);
    }

    /**
     * @see LangHelpers.applyPathEscapeAware()
     */
    private static Expression applyEscape(Expression node) {
        var cn = ClassHelper.makeWithoutCaching("nextflow.ast.LangHelpers");
        return callX(classX(cn), "applyPathEscapeAware", args(node));
    }

}
