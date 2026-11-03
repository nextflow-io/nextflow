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

import spock.lang.Specification

/**
 * Tests for the core-owned {@link AgentCallInfo} ThreadLocal seam (design §9.5/D6).
 */
class AgentCallInfoTest extends Specification {

    def cleanup() {
        AgentCallInfo.clear()
    }

    def 'should set and consume the resolved model round-trip'() {
        when:
        AgentCallInfo.setResolvedModel('gpt-4o-2024-08-06')
        then:
        AgentCallInfo.consumeResolvedModel() == 'gpt-4o-2024-08-06'
    }

    def 'consume should clear the value (second consume returns null)'() {
        given:
        AgentCallInfo.setResolvedModel('m1')
        expect:
        AgentCallInfo.consumeResolvedModel() == 'm1'
        AgentCallInfo.consumeResolvedModel() == null
    }

    def 'consume should be null-safe when nothing was set'() {
        expect:
        AgentCallInfo.consumeResolvedModel() == null
    }

    def 'should tolerate a null resolved model'() {
        when:
        AgentCallInfo.setResolvedModel(null)
        then:
        AgentCallInfo.consumeResolvedModel() == null
    }
}
