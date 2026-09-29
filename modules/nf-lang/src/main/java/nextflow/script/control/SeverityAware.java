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

import org.codehaus.groovy.control.messages.Message;
import org.codehaus.groovy.control.messages.SyntaxErrorMessage;
import org.codehaus.groovy.syntax.SyntaxException;

/**
 * Interface used by errors that should be reported as
 * warnings rather than failing the compilation.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public interface SeverityAware {
    boolean isSoftError();

    /**
     * Determine whether a compilation error is a soft error, i.e. it should
     * be reported as a warning instead of failing the compilation.
     *
     * @param message
     */
    static boolean isSoftError(Message message) {
        return message instanceof SyntaxErrorMessage sem && isSoftError(sem.getCause());
    }

    static boolean isSoftError(SyntaxException cause) {
        return cause instanceof SeverityAware sa && sa.isSoftError();
    }
}
