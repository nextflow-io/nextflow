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

import java.lang.reflect.Field

import groovyx.gpars.dataflow.DataflowReadChannel
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.runScript

/**
 * Regression test for the task-failure cascade bug.
 *
 * <p>When a tool's underlying process task hard-fails (exit ≠ 0) the session aborts the dataflow
 * network and interrupts the agent operator thread blocked on the tool output channel's
 * {@code .val} ({@link groovyx.gpars.dataflow.expression.DataflowExpression#getVal} throws a bare
 * {@link InterruptedException}). Before the fix, {@link ModuleToolBridge#call} swallowed that
 * {@link InterruptedException} into a {@code {"error":...}} tool result; langchain4j fed it back to
 * the model and the agent looped to {@code maxIterations}. The fix re-throws it as an
 * {@link AgentToolFatalError} (an {@link Error}, NOT an {@link Exception}) so it escapes
 * langchain4j's tool-execution {@code try/catch(Exception)} and aborts the run cleanly, while
 * restoring the thread's interrupt flag.
 *
 * <p>The test runs inside a live session (via {@code runScript} + a mock runner) so the bridge is
 * REAL and wired, then deterministically drives {@code call()} through the interrupt path (no
 * executor-timing dependency) by reflectively replacing the wired tool's captured output read
 * channel with a stub whose {@code getVal()} throws {@link InterruptedException} — exactly what a
 * session abort produces on the blocked operator thread.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(60)
class ModuleToolBridgeTaskFailureTest extends Dsl2Spec {

    def cleanup() {
        AgentRunnerProvider.testRunner = null
        Thread.interrupted() // clear any leaked interrupt flag on the test worker thread
    }

    /**
     * Reach into the bridge's private {@code tools} map and replace the named scalar tool's
     * captured output read channel ({@code toolOut}) with the given stub, so {@code call()}'s
     * blocking pull resolves to the stub.
     */
    private static void replaceToolOut(ModuleToolBridge bridge, String name, DataflowReadChannel stub) {
        final toolsField = ModuleToolBridge.getDeclaredField('tools')
        toolsField.setAccessible(true)
        final Map tools = (Map) toolsField.get(bridge)
        final tool = tools.get(name)
        assert tool != null
        final Field outField = tool.getClass().getDeclaredField('toolOut')
        outField.setAccessible(true)
        outField.set(tool, stub)
    }

    private static final String SCRIPT = '''
        nextflow.enable.types = true

        process greet {
            input:
            name: String

            output:
            greeting: String

            exec:
            greeting = "Hello ${name}!"
        }

        agent assistant {
            model 'm'
            instruction 'i'
            tools 'greet'

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

    def 'should abort (throw AgentToolFatalError) and restore interrupt flag when the blocking output pull is interrupted by task failure'() {
        given: 'a mock runner that, on the live operator thread, swaps the wired tool output for one that throws InterruptedException (== a session abort on task failure) and then invokes the tool'
        Throwable thrownByDispatch = null
        boolean interruptFlagAfter = false
        String dispatchResult = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            final bridge = (ModuleToolBridge) req.dispatch
            final throwing = Stub(DataflowReadChannel) {
                getVal() >> { throw new InterruptedException() }
            }
            replaceToolOut(bridge, 'greet', throwing)
            Thread.interrupted() // start from a clean interrupt status
            try {
                dispatchResult = req.dispatch.call('greet', '{"name":"Ada"}')
            }
            catch( Throwable t ) {
                thrownByDispatch = t
                interruptFlagAfter = Thread.currentThread().isInterrupted()
                Thread.interrupted() // clear so the operator can unwind cleanly
            }
            return 'done'
        } as AgentRunner

        when:
        runScript(SCRIPT)

        then: 'the dispatch did NOT return a {"error":...} result — it raised a fatal Error that aborts the run'
        dispatchResult == null
        thrownByDispatch instanceof AgentToolFatalError
        thrownByDispatch.message.contains('greet')
        thrownByDispatch.cause instanceof InterruptedException

        and: 'the thread interrupt flag was restored so the abort/teardown machinery observes it'
        interruptFlagAfter
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
