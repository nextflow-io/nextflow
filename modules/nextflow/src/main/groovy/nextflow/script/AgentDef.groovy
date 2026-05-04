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
import groovy.util.logging.Slf4j

/**
 * Runtime model for an `agent` definition.
 *
 * Mirrors {@link ProcessDef}'s contract — the script-level binding is callable
 * and chainable in a workflow — but execution is delegated to the future
 * `nf-agent` plugin runner. Calling {@link #run(Object[])} in this milestone
 * throws {@link UnsupportedOperationException} until the runner lands.
 *
 * @author Paolo Di Tommaso
 */
@Slf4j
@CompileStatic
class AgentDef extends BindableDef implements ChainableDef {

    static final String TYPE = 'agent'

    private BaseScript owner
    private String name
    private String simpleName
    private Closure body

    AgentDef(BaseScript owner, Closure body, String name) {
        this.owner = owner
        this.body = body
        this.name = name
        this.simpleName = name
    }

    @Override
    String getType() { TYPE }

    @Override
    String getName() { name }

    String getSimpleName() { simpleName }

    BaseScript getOwner() { owner }

    Closure getBody() { body }

    @Override
    ComponentDef cloneWithName(String name) {
        def copy = (AgentDef) this.clone()
        copy.@name = name
        copy.@simpleName = name.contains(':') ? name.tokenize(':').last() : name
        return copy
    }

    @Override
    Object run(Object[] args) {
        // Inherited path: BindableDef.invoke_a clones this and calls run().
        throw new UnsupportedOperationException('agent execution not yet implemented')
    }
}
