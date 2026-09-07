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
package nextflow.agent

import groovy.transform.CompileStatic

/**
 * Raised by the plugin's final structuring turn when the model refuses or returns
 * no usable content for the required structured-output schema (detected via a
 * {@code CONTENT_FILTER} finish reason or blank/null content).
 *
 * <p>A refusal is deterministic, so this fails fast with a clear, refusal-flavored
 * cause instead of a retry-storm. The exception is plugin-scoped; the core sees it
 * only as a {@link Throwable} through the {@code AgentRunner} SPI (boundary preserved).
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AgentStructuredOutputException extends RuntimeException {
    AgentStructuredOutputException(String message) {
        super(message)
    }
}
