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
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description
import nextflow.util.Duration

/**
 * Model the `agent` configuration scope. These settings are applied as
 * <em>defaults</em> for agents at runtime: an agent's own {@code model} /
 * {@code maxIterations} directives always take precedence; the scope only fills
 * in a value the agent did not declare.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ScopeName("agent")
@Description("""
    The `agent` scope controls the default settings applied to agents.
""")
@CompileStatic
class AgentConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        The default model (`provider/model`) used by an agent that does not declare a `model` directive.
    """)
    final String defaultModel

    @ConfigOption
    @Description("""
        The default maximum number of LLM iterations used by an agent that does not declare a `maxIterations` directive (default: `20`).
    """)
    final Integer maxIterationsDefault

    @ConfigOption
    @Description("""
        The amount of time to wait for an LLM chat request to complete before failing (default: `120 sec`).
    """)
    final Duration requestTimeout

    /* required by extension point -- do not remove */
    AgentConfig() {}

    AgentConfig(Map opts) {
        defaultModel = opts.defaultModel as String
        maxIterationsDefault = opts.maxIterationsDefault != null ? opts.maxIterationsDefault as Integer : null
        requestTimeout = opts.requestTimeout != null ? opts.requestTimeout as Duration : null
    }

    String getDefaultModel() { defaultModel }

    Integer getMaxIterationsDefault() { maxIterationsDefault }

    Duration getRequestTimeout() { requestTimeout }
}
