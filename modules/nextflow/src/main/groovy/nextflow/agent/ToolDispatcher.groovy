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
 * Callback the nf-agent plugin invokes to execute a tool: given the tool name
 * and the LLM-supplied arguments as a JSON string, it runs the corresponding
 * module/process and returns the tool result as a JSON string.
 *
 * Being a single-abstract-method interface, it can be supplied from Groovy as a
 * closure coerced {@code as ToolDispatcher}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
interface ToolDispatcher {
    String call(String toolName, String argsJson)
}
