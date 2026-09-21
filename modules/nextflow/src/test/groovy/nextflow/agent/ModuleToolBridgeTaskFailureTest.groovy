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

import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.runScript

/**
 * Regression test for the task-failure cascade bug.
 *
 * <p>When a tool's underlying process task hard-fails (exit ≠ 0) the session aborts the dataflow
 * network and interrupts the agent task thread blocked on its correlated reply variable.
 * Before the fix, {@link ModuleToolBridge#call} swallowed that
 * {@link InterruptedException} into a {@code {"error":...}} tool result; langchain4j fed it back to
 * the model and the agent looped to {@code maxIterations}. The fix re-throws it as an
 * {@link AgentToolFatalError} (an {@link Error}, NOT an {@link Exception}) so it escapes
 * langchain4j's tool-execution {@code try/catch(Exception)} and aborts the run cleanly, while
 * restoring the thread's interrupt flag.
 *
 * <p>The test runs a genuinely failing tool process inside a live session, exercising the
 * request queue, dynamic process invocation, session abort and reply interruption end to end.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(60)
class ModuleToolBridgeTaskFailureTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
        Thread.interrupted() // clear any leaked interrupt flag on the test worker thread
    }

    private static final String SCRIPT = '''
        nextflow.enable.types = true

        process greet {
            input:
            name: String

            output:
            greeting: String

            exec:
            throw new IllegalStateException('tool failed')
        }

        agent assistant {
            model 'm'
            instruction 'i'
            tools 'nf:module_run:greet'

            input:
            request: String

            output:
            answer: String

            prompt:
            """
            ${request}
            """
        }

        workflow {
            assistant(channel.of('hi')).view { it }
        }
        '''

    def 'should abort when the request-scoped tool process fails'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            req.dispatch.call('greet', '{"name":"Ada"}')
        } as AgentRunner

        when:
        Throwable failure = null
        try {
            runScript(SCRIPT)
        }
        catch( Throwable e ) {
            failure = e
        }

        then:
        failure != null
        failure.message.contains('tool failed') || failure.message.contains('greet')
    }

    def 'should still return a recoverable error result for a genuine dispatch-level failure (unknown tool)'() {
        given:
        String dispatchResult = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            dispatchResult = req.dispatch.call('nope', '{}')
            return 'done'
        } as AgentRunner

        when:
        runScript(SCRIPT)

        then: 'a dispatch-level error is recoverable: returned as a {"error":...} tool result, not thrown'
        dispatchResult != null
        new groovy.json.JsonSlurper().parseText(dispatchResult).error.contains('nope')
    }

}
