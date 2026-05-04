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
import nextflow.Session

/**
 * Factory for {@link AgentDef} instances. Counterpart to {@link ProcessFactory}.
 *
 * Deliberately thin in this milestone — the heavy lifting (chat-model construction,
 * tool registration, agent loop) lives in the future nf-agent plugin runner.
 *
 * @author Paolo Di Tommaso
 */
@Slf4j
@CompileStatic
class AgentFactory {

    private Session session
    private BaseScript owner

    AgentFactory(BaseScript ownerScript, Session session) {
        this.owner = ownerScript
        this.session = session
    }

    AgentDef newAgent(String name, Closure body) {
        new AgentDef(owner, body, name)
    }
}
