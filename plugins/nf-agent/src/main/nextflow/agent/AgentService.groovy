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
 * Minimal proxy interface implemented at runtime by langchain4j
 * {@code AiServices}. A single {@code chat} method drives one agent turn-set:
 * the proxy advertises the registered tools, dispatches tool-execution
 * requests, feeds results back to the model and returns the model's final
 * text answer.
 *
 * Named {@code AgentService} (not {@code Agent}) to avoid clashing with the
 * langchain4j-agentic {@code @Agent} annotation. The method carries no
 * annotations so the full user text is passed verbatim as the user message.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
interface AgentService {

    String chat(String userMessage)

}
