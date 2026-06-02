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

    static AgentRunner get() {
        if( testRunner != null )
            return testRunner
        final all = Plugins.getPriorityExtensions(AgentRunner)
        if( !all )
            throw new AbortOperationException("No agent runner available - enable the `nf-agent` plugin (e.g. add `plugins { id 'nf-agent' }` to your config)")
        return all.first()
    }
}
