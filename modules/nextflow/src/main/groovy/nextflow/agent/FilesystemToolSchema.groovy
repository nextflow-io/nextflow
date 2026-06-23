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
 * Descriptor for the generic {@code filesystem} agent tool. The work dir is bound
 * per-call from the dispatch context (never supplied by the LLM); the LLM provides
 * a path (relative to or inside the sandbox), an operation, and content for writes.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class FilesystemToolSchema {

    static final String NAME = 'filesystem'

    static ToolDescriptor descriptor() {
        final input = [
            type: 'object',
            properties: [
                path     : [type: 'string', description: 'File or directory path within the agent sandbox (the agent work dir) or a module-output path returned by a previous tool call.'],
                operation: [type: 'string', enum: ['read','write','list','exists'], description: 'The filesystem operation to perform.'],
                content  : [type: 'string', description: 'Content to write (required for the write operation; ignored otherwise).'],
            ],
            required: ['path','operation'],
            additionalProperties: false,
        ] as Map
        final desc = 'Read, write, list, or check files within the agent sandbox. Writes are confined to the agent work dir; reads may also target module-output paths returned by module_run. Small text files are returned inline; binary/large files are returned as path handles.'
        return new ToolDescriptor(NAME, desc, input, null)
    }
}
