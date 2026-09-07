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

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.agent.AgentLaunchSpec
import nextflow.agent.rpc.AgentRpcRegistration
import nextflow.agent.AgentRunner
import nextflow.agent.AgentRunnerRequest
import org.pf4j.Extension

/**
 * Pi SDK runner. Production execution is described by {@link #getLaunchSpec()} and runs as a
 * canonical Nextflow task in the runner container image, through the {@code agent-rpc} proxy.
 *
 * <p>The plugin ships no runtime of its own: both the proxy and the Node harness are provided by
 * the image, so an agent selecting this runner requires a container.
 */
@Extension
@CompileStatic
class PiAgentRunner implements AgentRunner {

    /**
     * Generated at build time from {@code build-image.sh} and this plugin's VERSION - see
     * {@code generateImageCoordinate} in build.gradle. Never hand-written, so the reference the
     * jar asks for is the reference the release publishes.
     */
    @PackageScope
    static final String IMAGE_RESOURCE = '/META-INF/nf-agent-pi-image.properties'

    private static final String IMAGE = parseImageCoordinate(PiAgentRunner.getResourceAsStream(IMAGE_RESOURCE))

    @Override
    String getName() { 'pi' }

    /**
     * The runner image, so an agent selecting {@code pi} needs no explicit {@code agent.container}.
     *
     * <p>Null when the generated resource is missing or unusable, which only a jar built outside
     * the Gradle build can be (PiAgentPackagingTest pins its presence in both artifacts). Null
     * rather than a throw on purpose: it degrades to core's existing "must declare a container"
     * failure, which is legible and recoverable by setting {@code agent.container}, instead of an
     * exception out of a static initialiser during PF4J's extension loading.
     */
    @Override
    String getDefaultContainer() { IMAGE }

    /** Read the {@code image} entry, or null if there is nothing usable to read. */
    @PackageScope
    static String parseImageCoordinate(InputStream stream) {
        if( stream == null )
            return null
        try {
            final props = new Properties()
            stream.withCloseable { InputStream it -> props.load(it) }
            return props.getProperty('image')?.trim() ?: null
        }
        catch( IOException e ) {
            return null
        }
    }

    /** The in-image command pair. Paths must match those the runner image is built with. */
    @Override
    AgentLaunchSpec getLaunchSpec() {
        return new AgentLaunchSpec(
            containerProxyCommand: ['/usr/local/bin/agent-rpc'],
            containerHarnessCommand: ['node', '/opt/nf-agent-pi/runner.mjs'] )
    }

    @Override
    AgentRpcRegistration register(AgentRunnerRequest request, boolean remote) {
        return AgentRpcBroker.get().register(request, remote)
    }

    /**
     * Not supported: this runner has no driver-local execution path. A canonical agent task runs
     * the container command described by {@link #getLaunchSpec()} and talks to the driver-side
     * broker over RPC, so nothing invokes the model from the driver JVM.
     */
    @Override
    String run(AgentRunnerRequest request) {
        throw new UnsupportedOperationException("Agent runner `${name}` runs as a container task - it cannot be invoked in the driver process")
    }
}
