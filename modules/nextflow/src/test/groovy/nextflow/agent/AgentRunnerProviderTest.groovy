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

import nextflow.exception.AbortOperationException
import spock.lang.Specification

class AgentRunnerProviderTest extends Specification {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    def 'should return the test runner when set'() {
        given:
        def runner = { AgentRunnerRequest req -> "echo:${req.prompt}" } as AgentRunner
        AgentRunnerProvider.testRunner = runner

        expect:
        AgentRunnerProvider.get().is(runner)
    }

    def 'should fail with a helpful error when no runner is available'() {
        when:
        AgentRunnerProvider.get()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('nf-agent')
    }
}
