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
package nextflow.script.dsl;

/**
 * DSL scope for agent definitions.
 *
 * Mirrors {@link ProcessDsl}: the outer interface is the definition scope
 * (where the agent's typed `input:` parameters are declared as locals), and
 * the nested {@link DirectiveDsl} declares the agent directive methods that
 * may appear at the top of an agent body.
 */
public interface AgentDsl extends DslScope {

    interface DirectiveDsl extends DslScope {

        @Description("""
            The `model` directive selects the LLM in `provider/model` form (e.g. `openai/gpt-5-mini`).
        """)
        void model(String value);

        @Description("""
            The `instruction` directive sets the agent system prompt (its role/persona).
        """)
        void instruction(String value);

        @Description("""
            The `tools` directive declares the modules the agent may invoke as tools.
        """)
        void tools(Object... values);

        @Description("""
            The `maxIterations` directive caps the LLM tool-calling loop.
        """)
        void maxIterations(Integer value);
    }
}
