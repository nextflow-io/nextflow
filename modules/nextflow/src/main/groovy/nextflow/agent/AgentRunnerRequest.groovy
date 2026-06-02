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

import groovy.transform.Canonical
import groovy.transform.CompileStatic

/**
 * Immutable request passed to an {@link AgentRunner}: the resolved model, the
 * system instruction, the rendered user prompt, the iteration cap, the
 * (currently unused) tool list for forward compatibility, the JSON schema
 * describing the expected structured output, and the input record serialized
 * as JSON.
 *
 * Being {@code @Canonical}, the positional constructor order is:
 * {@code (model, instruction, prompt, maxIterations, tools, outputSchema, inputJson)}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
@CompileStatic
class AgentRunnerRequest {
    String model
    String instruction
    String prompt
    int maxIterations
    List tools
    Map outputSchema
    String inputJson
}
