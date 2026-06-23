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
 * Fatal, non-recoverable error raised by {@link ModuleToolBridge#call} when the underlying
 * tool <i>process task</i> hard-fails (exit ≠ 0) and the session aborts the dataflow network,
 * interrupting the agent operator thread that is blocked on the tool's output channel.
 *
 * <p><b>Why an {@link Error} and not an {@link Exception}.</b> The agent dispatch runs inside a
 * langchain4j {@code AiServices} tool-execution loop. langchain4j wraps every {@code ToolExecutor}
 * call in a {@code try/catch(java.lang.Exception)} ({@code ToolService.executeWithErrorHandling}):
 * any thrown {@link Exception} — including a {@link RuntimeException} — is caught and converted
 * into an <i>error tool result</i> fed back to the model, which then keeps looping until the
 * iteration cap. To make a hard task failure abort the run cleanly (rather than loop to
 * {@code maxIterations}) the bridge must throw something OUTSIDE that {@code Exception} catch.
 * An {@link Error} propagates out of {@code agent.chat(...)} and reaches the agent operator's
 * {@code onException} listener, which calls {@code session.abort(t)} and poisons the tool queues.
 *
 * <p>This type deliberately carries NO langchain4j reference so the plugin boundary (no langchain4j
 * types under {@code modules/nextflow/src}) stays clean.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AgentToolFatalError extends Error {

    AgentToolFatalError(String message, Throwable cause) {
        super(message, cause)
    }

}
