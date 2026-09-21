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

import dev.langchain4j.agent.tool.ToolSpecification
import dev.langchain4j.model.chat.request.json.JsonObjectSchema
import groovy.transform.CompileStatic

/**
 * Maps a core, langchain4j-free {@link ToolDescriptor} onto a langchain4j
 * {@link ToolSpecification} so the LLM can be told which tools it may call.
 *
 * The descriptor's portable {@code inputSchema} {@link Map} (the same shape
 * {@link JsonSchemaMapper} consumes) becomes the tool's {@code parameters}
 * {@link JsonObjectSchema}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ModuleToolAdapter {

    /**
     * Build a langchain4j {@link ToolSpecification} from the given descriptor.
     *
     * @param d the portable tool descriptor (name, description, input schema)
     * @return a {@link ToolSpecification} whose {@code parameters} schema is
     *         derived from {@code d.inputSchema}
     */
    static ToolSpecification toToolSpecification(ToolDescriptor d) {
        if( !d )
            throw new IllegalArgumentException("Tool descriptor cannot be null")
        if( !d.name )
            throw new IllegalArgumentException("Tool descriptor `name` is required")
        final JsonObjectSchema params = JsonSchemaMapper.toObjectSchema(d.inputSchema ?: [type: 'object'])
        return ToolSpecification.builder()
            .name(d.name)
            .description(d.description)
            .parameters(params)
            .build()
    }
}
