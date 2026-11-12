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
import groovy.transform.PackageScope
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins

/**
 * Resolves the active {@link AgentRunner} from the loaded plugins. A package-scope
 * {@code testRunner} seam allows unit tests to inject a runner without a plugin.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AgentRunnerProvider {

    @PackageScope
    static AgentRunner testRunner

    @PackageScope
    static List<AgentRunner> testRunners

    /**
     * Resolve a runner by its stable extension name. When no name is supplied a
     * single installed runner is selected for backwards compatibility; multiple
     * runners are deliberately treated as ambiguous instead of depending on PF4J
     * extension ordering.
     */
    static AgentRunner get(String name = null) {
        if( testRunner != null )
            return testRunner
        final all = testRunners != null ? testRunners : Plugins.getPriorityExtensions(AgentRunner)
        if( !all )
            throw new AbortOperationException("No agent runner available - enable an agent runner plugin (for example `nf-agent` or `nf-agent-pi`)")
        if( name ) {
            final matches = all.findAll { it.getName() == name }
            if( matches.size() == 1 )
                return matches.first()
            if( matches.size() > 1 )
                throw new AbortOperationException("Multiple agent runner extensions use the name `${name}` - disable the duplicate plugin")
            throw new AbortOperationException("Unknown agent runner `${name}` - available runners: ${availableNames(all)}")
        }
        if( all.size() == 1 )
            return all.first()
        throw new AbortOperationException("Multiple agent runners are available (${availableNames(all)}) - select one with `agent.runner` in nextflow.config")
    }

    private static String availableNames(List<AgentRunner> all) {
        all.collect { it.getName() }.toSorted().join(', ')
    }
}
