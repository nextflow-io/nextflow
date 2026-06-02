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
 * Portable, langchain4j-free descriptor of a module/process exposed to the LLM
 * as a tool. The input and output schemas are plain JSON-schema {@link Map}s
 * (same shape produced by {@link RecordSchema} and {@link ProcessToolSchema}) so
 * that core stays free of any LLM client dependency; the nf-agent plugin maps
 * them onto langchain4j's {@code ToolSpecification}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
@CompileStatic
class ToolDescriptor {
    String name
    String description
    Map inputSchema
    Map outputSchema
}
