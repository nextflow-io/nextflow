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
package nextflow.trace

import groovy.transform.CompileStatic

/**
 * Writes one terminal activity report.
 *
 * The bytes belong to the implementation. {@link Osc7501Emitter} is the
 * Program Status Protocol (OSC 7501) revision 0.4. A later revision of that
 * proposal can replace this without touching the run-state tracker.
 */
@CompileStatic
interface TerminalActivityEmitter {

    /**
     * @param state protocol state token, such as {@code idle} or {@code working}
     * @param message one human-readable line, or {@code null}
     */
    void emit(String state, String message)

}
