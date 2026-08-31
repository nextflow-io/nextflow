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

import groovy.transform.NamedParam;

/**
 * DSL scope for agent definitions.
 *
 * Mirrors {@link ProcessDsl}: the outer interface is the definition scope
 * (where the agent's typed `input:` parameters are declared as locals), the
 * nested {@link DirectiveDsl} declares the agent directive methods that
 * may appear at the top of an agent body, and {@link AgentOutputDsl} the functions
 * available in the `output:` section.
 */
public interface AgentDsl extends DslScope {

    interface DirectiveDsl extends DslScope {

        @Description("""
            The `goal` directive states a high-level objective that steers the agent's multi-turn loop. It is advisory: the model is encouraged to keep working until the goal is met, while `maxIterations` remains the hard cap.
        """)
        void goal(String value);

        @Description("""
            The `instruction` directive sets the agent system prompt (its role/persona).
        """)
        void instruction(String value);

        @Description("""
            The `label` directive annotates the agent with a mnemonic identifier, which can be used to apply configuration in the `agent` scope through a `withLabel:` selector. It can be specified more than once.
        """)
        void label(String value);

        @Description("""
            The `maxIterations` directive caps the LLM tool-calling loop.
        """)
        void maxIterations(Integer value);

        @Description("""
            The `model` directive selects the LLM in `provider/model` form (e.g. `openai/gpt-5-mini`).
        """)
        void model(String value);

        @Description("""
            The `skills` directive declares the agent skills (SKILL.md folders) the agent may use. Each entry is a local skill name (resolved under the `skills/` directory beside the script) or a remote GitHub reference (`github.com/<org>/<repo>[@rev]`) cloned and cached into that same `skills/` directory.
        """)
        void skills(Object... values);

        @Description("""
            The `tools` directive declares the tools the agent may invoke, as namespaced references of the form `family[:group]:name`. `nf:module_run` exposes every in-scope module or process, `nf:module_run:<NAME>` a single one and `nf:module_run:<PREFIX>*` those matching; `fs:*` selects the six sandboxed file tools (`read`, `write`, `edit`, `ls`, `grep`, `find`); `shell:bash` a shell inside the runner container (`pi` runner only). An agent declaring no tools gets none, and a reference matching nothing is an error.
        """)
        void tools(Object... values);
    }

    /**
     * Functions available in an agent `output:` section.
     *
     * Deliberately a SUBSET of {@link ProcessDsl.OutputDslV2}, and a subset BY CONSTRUCTION: both
     * scopes inherit `file`/`files` from the shared {@link FileOutputDsl}, so an option added to
     * one is added to the other. `eval`/`stdout` stay process-only because an agent has no task
     * script to read them back from — leaving them undeclared makes them a resolution error
     * instead of a runtime surprise.
     *
     * <p>Named {@code AgentOutputDsl} rather than {@code OutputDsl} to stay distinct from the
     * unrelated top-level {@link OutputDsl} (the workflow output DSL) in this same package.
     */
    /** Functions available in an agent `output:` section; `file`/`files` are shared with the process scope. */
    interface AgentOutputDsl extends FileOutputDsl {
    }
}
