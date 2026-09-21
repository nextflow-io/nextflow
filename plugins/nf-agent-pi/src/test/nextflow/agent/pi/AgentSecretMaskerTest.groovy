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
package nextflow.agent.pi

import nextflow.SysEnv
import spock.lang.Specification

/**
 * Pins design D9: redaction keys off the RESOLVED credential, not off an environment variable
 * name. A key from `agent.apiKey` (typically `secrets.LLM_KEY`) has no variable name, so the
 * name-driven sweep that predates this cannot find it.
 */
class AgentSecretMaskerTest extends Specification {

    def 'redacts a credential that has NO environment variable name'() {
        given: 'the D9 case: the key came from `agent.apiKey`, so an env sweep cannot discover it'
        SysEnv.push([:])
        and: 'and it is NOT OpenAI-shaped, so the `sk-` backstop pattern cannot catch it either'
        def key = 'gateway-credential-4a71b2'

        expect:
        AgentSecretMasker.redact("401 unauthorized for key ${key}".toString(), key) ==
            '401 unauthorized for key [REDACTED]'

        and: 'every occurrence goes, not just the first'
        AgentSecretMasker.redact("${key} and again ${key}".toString(), key) == '[REDACTED] and again [REDACTED]'

        and: 'without the resolved key the same text leaks -- this is the regression being fixed'
        AgentSecretMasker.redact("401 unauthorized for key ${key}".toString()) ==
            "401 unauthorized for key ${key}"

        cleanup:
        SysEnv.pop()
    }

    def 'sweeps the exported provider variables through SysEnv, not System.getenv'() {
        given: 'the sweep is kept as a backstop for a key the user still exports'
        SysEnv.push([OPENAI_API_KEY: 'sk-exported-01234567', ANTHROPIC_API_KEY: 'ant-exported-9999', NXF_AGENT_API_KEY: 'nxf-exported-8888'])

        expect: 'reading through SysEnv is what lets a test swap the environment at all'
        AgentSecretMasker.redact('failed with ant-exported-9999') == 'failed with [REDACTED]'
        AgentSecretMasker.redact('failed with nxf-exported-8888') == 'failed with [REDACTED]'

        and: 'an unrelated variable is left alone'
        AgentSecretMasker.redact('failed at http://localhost:8000/v1') == 'failed at http://localhost:8000/v1'

        cleanup:
        SysEnv.pop()
    }

    def 'masks bearer headers and credential-shaped tokens the provider echoes back'() {
        given: 'the last resort: a credential that is neither resolved nor exported here'
        SysEnv.push([:])

        expect:
        AgentSecretMasker.redact('Authorization: Bearer abcd1234efgh') == 'Authorization: Bearer [REDACTED]'
        AgentSecretMasker.redact('rejected sk-abcdef0123456789') == 'rejected [REDACTED]'
        AgentSecretMasker.redact('rejected rk-abcdef0123456789') == 'rejected [REDACTED]'
        AgentSecretMasker.redact('rejected pk-abcdef0123456789') == 'rejected [REDACTED]'

        and: 'a short `sk-` fragment is not a key and is left readable'
        AgentSecretMasker.redact('sk-short') == 'sk-short'

        cleanup:
        SysEnv.pop()
    }

    def 'is null-safe and leaves clean text untouched'() {
        given:
        SysEnv.push([:])

        expect: 'null in, null out -- callers pass optional message fields straight through'
        AgentSecretMasker.redact(null) == null
        AgentSecretMasker.redact(null, 'sk-whatever') == null
        and: 'an empty resolved key must not turn every character boundary into [REDACTED]'
        AgentSecretMasker.redact('provider refused request', '') == 'provider refused request'
        AgentSecretMasker.redact('provider refused request', null) == 'provider refused request'

        cleanup:
        SysEnv.pop()
    }
}
