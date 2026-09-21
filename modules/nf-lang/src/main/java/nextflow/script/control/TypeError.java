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
package nextflow.script.control;

import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.syntax.SyntaxException;

/**
 * Error reported by the type checking phase. A soft error is reported as a
 * warning rather than failing the compilation.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class TypeError extends SyntaxException implements PhaseAware, SeverityAware {

    private final boolean softError;

    public TypeError(String message, ASTNode node, boolean softError) {
        super(message, node);
        this.softError = softError;
    }

    public TypeError(String message, ASTNode node) {
        this(message, node, false);
    }

    @Override
    public int getPhase() {
        return Phases.TYPE_CHECKING;
    }

    @Override
    public boolean isSoftError() {
        return softError;
    }

}
