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

package nextflow.agent;

import java.lang.reflect.Proxy;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;

import nextflow.Global;
import nextflow.NF;
import nextflow.Session;
import nextflow.processor.TaskProcessor;
import nextflow.script.ScriptBinding;
import nextflow.script.ScriptLoaderFactory;
import nextflow.script.ScriptMeta;
import org.codehaus.groovy.runtime.InvokerHelper;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;
import org.junit.jupiter.api.io.TempDir;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

/**
 * Checks how the task of an agent tool call completes the call: a failed task that is
 * not ignored is fatal, and every path of the task completes the call.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Timeout(value = 60, unit = TimeUnit.SECONDS)
class ModuleToolTaskCompletionTest {

    @TempDir
    Path workDir;

    private final AtomicReference<Object> outcome = new AtomicReference<>();

    @BeforeEach
    void setup() {
        TaskProcessor.reset();
        ScriptMeta.reset();
        Global.reset();
        NF.init();
        // record the result or the error of a tool call, and fail the agent task with the error
        // a proxy, as the default methods of the Groovy interface are abstract for Java
        AgentRunnerProvider.testRunner = (AgentRunner) Proxy.newProxyInstance(getClass().getClassLoader(), new Class<?>[] { AgentRunner.class }, (proxy, method, args) -> {
            switch( method.getName() ) {
                case "run": break;
                case "getName": case "toString": return "test";
                case "hashCode": return System.identityHashCode(proxy);
                case "equals": return proxy == args[0];
                default: return null;
            }
            try {
                final var result = ((AgentRunnerRequest) args[0]).getDispatch().call("greet", "{\"name\":\"Ada\"}");
                outcome.set(result);
                return result;
            }
            catch( Throwable e ) {
                outcome.set(e);
                throw e;
            }
        });
    }

    @AfterEach
    void cleanup() {
        AgentRunnerProvider.testRunner = null;
        Thread.interrupted();
    }

    /**
     * Run a script on the local executor, whose native task handler reports the error of a
     * task like a real executor
     */
    private void run(String text) throws Throwable {
        final var session = new Session(new HashMap<>(Map.of("workDir", workDir.toString())));
        session.setBinding(new ScriptBinding());
        session.init(null, null, null, null);
        session.start();
        try {
            final var loader = ScriptLoaderFactory.create(session);
            // parse(String) is not part of the ScriptLoader interface
            InvokerHelper.invokeMethod(loader, "parse", text);
            loader.runScript();
            session.fireDataflowNetwork(false);
            session.await();
        }
        finally {
            session.destroy();
        }
        if( session.getError() != null )
            throw session.getError();
    }

    private static String script(String greet) {
        return """
            nextflow.enable.types = true

            process greet {
                %s
            }

            agent assistant {
                model 'm'
                instruction 'i'
                tools 'nf:module_run:greet'

                input:
                request: String

                output:
                stdout()

                prompt:
                \"""
                ${request}
                \"""
            }

            workflow {
                assistant(channel.of('hi')).view { it }
            }
            """.formatted(greet);
    }

    @Test
    void shouldFailTheToolCallWhenTheTaskFails() throws Throwable {
        final var text = script("""
            input:
            name: String

            output:
            greeting: String = 'unreachable'

            exec:
            throw new IllegalStateException('tool failed')
            """);

        assertThatThrownBy(() -> run(text));
        assertThat(outcome.get()).isInstanceOf(AgentToolFatalError.class);
    }

    @Test
    void shouldReturnAnErrorResultWhenTheTaskErrorIsIgnored() throws Throwable {
        final var text = script("""
            errorStrategy 'ignore'

            input:
            name: String

            output:
            greeting: String = 'unreachable'

            exec:
            throw new IllegalStateException('tool failed')
            """);

        run(text);
        assertThat(outcome.get()).asString().contains("\"error\"").contains("ignored");
    }

    @Test
    void shouldFailTheToolCallWhenTheWhenGuardFails() throws Throwable {
        final var text = script("""
            input:
            name: String

            output:
            greeting: String = "hello ${name}"

            when:
            name.foo()

            exec:
            true
            """);

        assertThatThrownBy(() -> run(text));
        assertThat(outcome.get()).isInstanceOf(AgentToolFatalError.class);
    }

    @Test
    void shouldReturnTheOutputsOfTheTask() throws Throwable {
        final var text = script("""
            input:
            name: String

            output:
            greeting: String = "hello ${name}"

            exec:
            true
            """);

        run(text);
        assertThat(outcome.get()).asString().contains("hello Ada");
    }
}
