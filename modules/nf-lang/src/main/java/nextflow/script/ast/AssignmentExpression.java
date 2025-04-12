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

import org.codehaus.groovy.ast.expr.BinaryExpression;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.syntax.Token;

import static org.codehaus.groovy.ast.tools.GeneralUtils.ASSIGN;

/**
 * Convenience class for assignment expressions.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class AssignmentExpression extends BinaryExpression {

    public AssignmentExpression(Expression target, Token op, Expression value) {
        super(target, op, value);
    }

    public AssignmentExpression(Expression target, Expression value) {
        super(target, ASSIGN, value);
    }
}
