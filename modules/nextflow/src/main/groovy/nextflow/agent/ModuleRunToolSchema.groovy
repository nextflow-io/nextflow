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
 * Builds the single {@code module_run} tool descriptor. The LLM picks a module by
 * name (constrained to an enum of the script's included modules) and supplies a
 * generic {@code args} object marshalled per that module's existing binding rules.
 * Per-module input hints are aggregated into the description so the LLM knows what
 * each module needs without a separate tool per module.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ModuleRunToolSchema {

    static final String NAME = 'module_run'

    static ToolDescriptor descriptor(List<String> moduleNames, Map<String,String> hints) {
        final input = [
            type: 'object',
            properties: [
                module: [type: 'string', enum: new ArrayList<String>(moduleNames),
                         description: 'Which module to run. Must be one of the available modules.'],
                args  : [type: 'object', additionalProperties: true,
                         description: 'The module inputs, as a flat object matching the chosen module (see the per-module input hints in this tool description).'],
            ],
            required: ['module','args'],
            additionalProperties: false,
        ] as Map
        final sb = new StringBuilder()
        sb.append('Run one of the available Nextflow modules as a real task and return its outputs as JSON. Available modules and their inputs:')
        for( final name : moduleNames ) {
            sb.append('\n\n### ').append(name).append('\n')
            sb.append(hints.get(name) ?: 'inputs: (see module documentation)')
        }
        sb.append('\n\nFile/path outputs are returned as absolute path strings (read them with the filesystem tool if needed), never file contents.')
        return new ToolDescriptor(NAME, sb.toString(), input, null)
    }
}
