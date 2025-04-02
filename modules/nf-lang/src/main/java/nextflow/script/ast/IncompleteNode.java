/*
 * Copyright 2024-2025, Seqera Labs
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
package nextflow.script.ast;

import org.codehaus.groovy.ast.expr.EmptyExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;

/**
 * An incomplete script declaration, used to provide more
 * contextual error messages and completions.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class IncompleteNode extends ExpressionStatement {
    public final String text;

    public IncompleteNode(String text) {
        super(EmptyExpression.INSTANCE);
        this.text = text;
    }

    public IncompleteNode(Expression expression) {
        super(expression);
        this.text = null;
    }
}
