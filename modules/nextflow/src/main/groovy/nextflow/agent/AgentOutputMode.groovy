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
 * How an agent's answer comes back: as the model's plain text, as a scalar wrapped in the
 * declared output name, as a record, or as an object holding one entry per declared output.
 *
 * <p>The single fact {@link AgentOutputPlan} is built around -- both the decoding of a canonical
 * task's terminal frame and the binding of an in-JVM runner's result are functions of it.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
enum AgentOutputMode {
    TEXT,
    SCALAR_CONTRACT,
    RECORD,
    WRAPPED
}
