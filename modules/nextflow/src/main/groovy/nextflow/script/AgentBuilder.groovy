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
package nextflow.script

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j

/**
 * Runtime delegate that captures an agent body when the lowered agent closure
 * runs. Directives are captured via {@link #methodMissing} against a known set;
 * inputs/outputs via {@code _input_}/{@code _output_}; the prompt via
 * {@link #withPrompt}. {@link #build} produces the populated {@link AgentDef}.
 *
 * Distinct from {@code nextflow.script.dsl.AgentDsl} (the compile-time
 * resolution scope).
 */
@Slf4j
@CompileStatic
class AgentBuilder {

    static final List<String> DIRECTIVES = ['model', 'instruction', 'goal', 'tools', 'maxIterations']

    private BaseScript ownerScript
    private String agentName

    private final Map<String,Object> directives = new LinkedHashMap<>()
    private final List<AgentInput> inputs = new ArrayList<>()
    private final List<AgentOutput> outputs = new ArrayList<>()
    private PromptDef prompt

    AgentBuilder(BaseScript ownerScript, String agentName) {
        this.ownerScript = ownerScript
        this.agentName = agentName
    }

    protected void checkName(String name) {
        if( !DIRECTIVES.contains(name) )
            throw new IllegalArgumentException("Unknown agent directive `${name}`")
    }

    @PackageScope
    Object methodMissing(String name, Object args) {
        checkName(name)
        final values = args instanceof Object[] ? (args as List) : [args]
        directives.put(name, values.size() == 1 ? values[0] : values)
        return null
    }

    void _input_(String name, Class type) {
        inputs.add(new AgentInput(name, type))
    }

    void _output_(String name, Class type) {
        outputs.add(new AgentOutput(name, type))
    }

    AgentBuilder withPrompt(PromptDef prompt) {
        this.prompt = prompt
        return this
    }

    AgentDef build() {
        if( prompt == null )
            throw new IllegalStateException("Missing prompt in agent `${agentName}` definition")
        return new AgentDef(ownerScript, agentName, directives, inputs, outputs, prompt)
    }

    @CompileStatic
    static class AgentInput {
        final String name
        final Class type
        AgentInput(String name, Class type) { this.name = name; this.type = type }
    }

    @CompileStatic
    static class AgentOutput {
        final String name
        final Class type
        AgentOutput(String name, Class type) { this.name = name; this.type = type }
    }
}
