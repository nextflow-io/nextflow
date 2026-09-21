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

import java.util.zip.ZipFile

import spock.lang.Specification

/**
 * Pins the plugin's DISTRIBUTION shape: the agent proxy and the Node harness are shipped in the
 * runner container image, so neither the jar nor the distribution zip carries them. Re-vendoring
 * them -- the 186 MB, host-arch-only artifact this replaced -- fails here rather than at release
 * time.
 *
 * <p>Also pins the seam between the two halves of that split: the in-image paths
 * {@link PiAgentRunner#getLaunchSpec} emits are only correct because the Dockerfile puts the proxy
 * and the harness exactly there. Nothing else couples the two files, and a drift is a `No such
 * file` inside the container on every agent task.
 */
class PiAgentPackagingTest extends Specification {

    /**
     * The generated image coordinate as an ARCHIVE ENTRY name - {@link PiAgentRunner#IMAGE_RESOURCE}
     * names it the way a classpath lookup does, with a leading slash. Derived from that constant so
     * renaming the resource cannot leave this spec asserting on a name nothing reads.
     */
    private static final String IMAGE_RESOURCE = PiAgentRunner.IMAGE_RESOURCE.substring(1)

    /** An artifact Gradle just built; its path is handed to the test task (see build.gradle). */
    private static File artifact(String property) {
        final path = System.getProperty(property)
        assert path, "${property} is unset - run this spec through Gradle (:plugins:nf-agent-pi:test)"
        final file = new File(path)
        assert file.exists(), "the artifact was not built: ${file}"
        return file
    }

    private static List<String> entriesOf(File archive) {
        return new ZipFile(archive).withCloseable { ZipFile zip -> zip.entries().toList().collect { it.name } }
    }

    def 'the plugin jar ships no vendored runtime'() {
        given:
        final jar = artifact('nf.agent.pi.jar')

        when:
        final entries = entriesOf(jar)

        then: 'it IS the plugin jar - asserted first, because an empty jar would satisfy everything below'
        entries.contains('nextflow/agent/pi/PiAgentRunner.class')
        entries.contains('nextflow/agent/pi/AgentRpcBroker.class')

        and: 'nothing under the prefix the deleted PiRuntime used to extract to a temp dir'
        entries.every { !it.startsWith('pi-runtime/') }

        and: 'each of these was a defect of its own: the 173 MB Node tree, the arch-specific Go proxy, the harness'
        entries.every { !it.contains('node_modules') }
        entries.every { !(it == 'agent-rpc' || it.endsWith('/agent-rpc')) }
        entries.every { !it.endsWith('runner.mjs') }

        and: 'the generated image coordinate DOES ship - without it the runner declares no container'
        entries.contains(IMAGE_RESOURCE)

        and: 'so the jar is Groovy classes plus the plugin metadata - measured ~37 KB'
        jar.length() < 128 * 1024
    }

    def 'the distribution zip ships no vendored runtime either'() {
        given: 'the artifact a user actually downloads, not just the jar inside it'
        final zip = artifact('nf.agent.pi.zip')

        when:
        final entries = entriesOf(zip)

        then: 'it IS the plugin distribution'
        entries.any { it.endsWith('nextflow/agent/pi/PiAgentRunner.class') }

        and: 'the generated image coordinate reached the artifact a user installs, not just the jar'
        // asserted here as well as on the jar for the reason this file already gives in reverse:
        // the two artifacts are assembled by different tasks, so a resource can reach one alone
        entries.any { it.endsWith(IMAGE_RESOURCE) }

        and: 'no runtime re-vendored through packagePlugin, which a jar-only assertion would miss'
        entries.every { !it.contains('node_modules') }
        entries.every { !it.contains('pi-runtime/') }
        entries.every { !it.endsWith('/agent-rpc') }
        entries.every { !it.endsWith('runner.mjs') }

        and: 'BouncyCastle actually SHIPS -- the broker builds its per-run TLS identity with it, and weakening `api` to `compileOnly` in build.gradle would pass every other assertion here and then fail as NoClassDefFoundError on the first agent task of every run'
        libJars(entries).containsAll(['bcprov-jdk18on', 'bcpkix-jdk18on'])

        and: 'and the runtime is pinned by IDENTITY, not by weight'
        libJars(entries) == EXPECTED_LIB_JARS
    }

    /**
     * The dependency jars the distribution is expected to carry, versions stripped.
     *
     * <p>This replaces a byte ceiling, which had quietly stopped guarding anything: it was written
     * when the only weight was the gRPC transport, then BouncyCastle added ~10.8 MB and left ~1.5 MB
     * of headroom under a 24 MB bound. A ceiling degrades silently as things slip under it, and it
     * answers the wrong question -- a re-vendored runtime is a defect at ANY size, while a bcprov
     * patch bump is harmless and would have failed the build with a message reading like a
     * re-vendoring regression.
     *
     * <p>Identity does not degrade. Versions are stripped so an ordinary bump is invisible here;
     * only a jar APPEARING or DISAPPEARING fails, which is exactly the change that deserves a human
     * decision. If you added a dependency on purpose, add it here in the same commit.
     */
    private static final Set<String> EXPECTED_LIB_JARS = [
        // the broker's TLS identity -- see AgentRpcTlsCredentials
        'bcpkix-jdk18on', 'bcprov-jdk18on', 'bcutil-jdk18on',
        // the agent RPC transport
        'grpc-api', 'grpc-context', 'grpc-core', 'grpc-netty-shaded', 'grpc-stub', 'grpc-util',
        'perfmark-api', 'gson', 'guava', 'failureaccess', 'listenablefuture',
        // annotation-only artifacts the above drag in
        'animal-sniffer-annotations', 'annotations', 'checker-qual', 'error_prone_annotations',
        'j2objc-annotations', 'jsr305',
    ] as Set<String>

    /**
     * Artifact names of the jars under {@code lib/}, with the trailing {@code -<version>} removed.
     * The pattern anchors on the first hyphen followed by a DIGIT, so `grpc-netty-shaded-1.75.0.jar`
     * and `listenablefuture-9999.0-empty-to-avoid-conflict-with-guava.jar` both reduce to their
     * artifact name rather than being truncated at the first hyphen.
     */
    private static Set<String> libJars(List<String> entries) {
        return entries
            .findAll { it.startsWith('lib/') && it.endsWith('.jar') }
            .collect { it.substring('lib/'.length()).replaceFirst(/-\d[^\/]*\.jar$/, '') } as Set<String>
    }

    def 'the launch spec paths are the paths the runner image is built with'() {
        given:
        final spec = new PiAgentRunner().getLaunchSpec()
        final dockerfile = artifact('nf.agent.pi.dockerfile').readLines()

        when: 'the destinations the image build puts the proxy and the harness at'
        final workDir = dockerfile.findAll { it.trim().startsWith('WORKDIR ') }.last().trim() - 'WORKDIR '
        final copied = dockerfile
            .findAll { it.trim().startsWith('COPY ') }
            .collectEntries { line ->
                final args = line.trim().split(/\s+/).findAll { !it.startsWith('--') && it != 'COPY' }
                [ (args[-2]): resolveIn(workDir, args[-1]) ]
            }

        then: 'the proxy binary the launch command execs'
        spec.containerProxyCommand == [ copied['/agent-rpc'] ]

        and: 'the harness the proxy runs after the `--` separator'
        spec.containerHarnessCommand == [ 'node', copied['harness/runner.mjs'] ]
    }

    def 'the declared image is the image build-image.sh publishes'() {
        given:
        final declared = new PiAgentRunner().getDefaultContainer()
        final script = artifact('nf.agent.pi.buildscript')

        when: 'ask the publishing script itself, which needs no docker for `ref`'
        // NF_AGENT_PI_REGISTRY is a supported override (build-image.sh), and Groovy's no-arg
        // .execute() inherits the JVM environment, so a developer or a CI runner with it exported
        // would get a red test claiming the jar's coordinate is wrong. Hand the child an
        // environment with that one variable removed: what the jar embeds is DEFAULT_REGISTRY.
        final env = System.getenv()
                          .findAll { k, v -> k != 'NF_AGENT_PI_REGISTRY' }
                          .collect { k, v -> "$k=$v" as String }
        final proc = [ 'bash', script.absolutePath, 'ref' ].execute(env, null)
        final out = new StringBuffer(), err = new StringBuffer()
        // drained rather than read through proc.text: the script writes diagnostics to stderr,
        // and an unread stderr pipe is a deadlock waiting for a longer message
        proc.waitForProcessOutput(out, err)

        then: 'the two compositions - Gradle reading the script, and the script itself - agree'
        proc.exitValue() == 0
        declared
        declared == out.toString().trim()
    }

    /** Resolve a Dockerfile COPY destination, which may be relative to the current WORKDIR. */
    private static String resolveIn(String workDir, String destination) {
        if( destination.startsWith('/') )
            return destination
        final relative = destination.startsWith('./') ? destination.substring(2) : destination
        return "${workDir}/${relative}".toString()
    }
}
