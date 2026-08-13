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

import nextflow.processor.TaskConfig
import spock.lang.Specification

/**
 * The capability token must not leave the work directory. It is a bearer credential FOR the
 * provider API key ever since the driver started answering a `connect` with the key on the start
 * frame, and {@code TraceRecord.script} is persisted in the resume cache and POSTed to Seqera
 * Platform.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AgentTaskScriptTest extends Specification {

    /** The exact shape AgentDef.createCanonicalBody produces: every argument shell-quoted. */
    private static String launchScript(String token) {
        return "exec '/opt/nf-agent/agent-rpc' '--endpoint' 'driver.internal:41235' " +
            "'--invocation' 'inv-7f3a' '--fingerprint' 'a1b2c3' '--token' '${token}' " +
            "'--' '/usr/bin/node' '/opt/nf-agent/runner.mjs'"
    }

    private static AgentTaskInfo agentInfo() {
        return new AgentTaskInfo('pi', 'openai/gpt-5-mini', 'be helpful', null, 'do ${x}', 20, null, null, null)
    }

    def 'the capability token is redacted and nothing else is'() {
        given:
        def token = 'Zm9vYmFyLXNlY3JldC10b2tlbi12YWx1ZQ'

        when:
        def redacted = AgentTaskScript.redactCapabilityToken(launchScript(token))

        then: 'the token value is gone, in every spelling of itself'
        !redacted.contains(token)
        redacted.contains("'--token' '${AgentTaskScript.REDACTED}'")

        and: 'everything a recorded script is worth having for survives'
        // a fingerprint is a PUBLIC commitment, not a secret -- see docs/agent.mdx
        redacted.contains("'--endpoint' 'driver.internal:41235'")
        redacted.contains("'--invocation' 'inv-7f3a'")
        redacted.contains("'--fingerprint' 'a1b2c3'")
        redacted.contains('/opt/nf-agent/runner.mjs')
    }

    def 'an unquoted spelling is redacted too, as a backstop against a future call site'() {
        expect:
        AgentTaskScript.redactCapabilityToken('exec proxy --token abc123 --') ==
            "exec proxy --token ${AgentTaskScript.REDACTED} --"
        AgentTaskScript.redactCapabilityToken('exec proxy --token=abc123') ==
            "exec proxy --token=${AgentTaskScript.REDACTED}"
    }

    def 'a script with no token is returned unchanged, not reformatted'() {
        expect:
        AgentTaskScript.redactCapabilityToken(SCRIPT) === SCRIPT

        where:
        SCRIPT << [
            'echo hello world',
            "exec '/opt/nf-agent/agent-rpc' '--endpoint' 'x' '--insecure'",
            '',
            null ]
    }

    def 'only an agent task is redacted; every other task is byte-identical'() {
        given: 'the guard that keeps this off the path of every ordinary process task'
        def script = launchScript('a-token-value')

        expect: 'a process task passes through untouched -- same instance, not merely equal'
        AgentTaskScript.forTrace(new TaskConfig([tag: 'x']), script) === script
        AgentTaskScript.forTrace(null, script) === script

        and: 'a FORGED agentInfo directive is inert: the guard is instanceof, like LinObserver'
        AgentTaskScript.forTrace(new TaskConfig([agentInfo: 'not-an-agent']), script) === script

        and: 'while a real agent task is redacted'
        !AgentTaskScript.forTrace(new TaskConfig([agentInfo: agentInfo()]), script).contains('a-token-value')
    }

    def 'isAgentTask recognizes the carrier the lineage observer uses'() {
        expect:
        AgentTaskScript.isAgentTask(new TaskConfig([(AgentTaskInfo.CONFIG_KEY): agentInfo()]))
        !AgentTaskScript.isAgentTask(new TaskConfig([:]))
        !AgentTaskScript.isAgentTask(null)
    }
}
