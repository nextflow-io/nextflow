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

import groovy.transform.Canonical
import groovy.transform.CompileStatic

/**
 * Portable command description for a runner that supports canonical executor tasks.
 *
 * <p>A runner describes only HOW to launch itself, with the ABSOLUTE paths its proxy and harness
 * have INSIDE its runner image: a canonical agent task always runs in a container, so there is a
 * single command pair and no driver-local variant. The image itself is NOT part of this spec: it
 * comes from the `agent.container` directive when set -- exactly as `process.container` does for a
 * process -- and otherwise from {@link AgentRunner#getDefaultContainer}, which lets a runner that
 * ships its runtime in an image name that image itself.
 * {@link nextflow.agent.AgentLaunchConditions#requireContainerized} still rejects a configuration that would
 * not containerize the task, whichever of the two supplied the image, so these paths are never
 * used on the driver.
 */
@Canonical
@CompileStatic
class AgentLaunchSpec {
    List<String> containerProxyCommand
    List<String> containerHarnessCommand

    /**
     * Compose a launch command with arguments passed to the proxy before the proxy/harness
     * separator. Keeping the composition here avoids callers having to build and then parse the
     * command to find the separator.
     */
    List<String> command(List<String> proxyArgs = Collections.emptyList()) {
        if( !containerProxyCommand || !containerHarnessCommand )
            throw new IllegalStateException("Agent runner does not provide a container launch command")
        final result = new ArrayList<String>(containerProxyCommand)
        result.addAll(proxyArgs)
        result.add('--')
        result.addAll(containerHarnessCommand)
        return result
    }

    /**
     * The same command as a single POSIX shell line -- what a canonical agent task's script IS.
     * Kept here so the object that composes the command also owns how it is spelled, rather than
     * leaving each caller to quote the arguments it produced.
     */
    String shellCommand(List<String> proxyArgs) {
        return 'exec ' + command(proxyArgs).collect { quote(it) }.join(' ')
    }

    /** POSIX shell quoting for generated canonical agent commands. */
    private static String quote(Object value) {
        final String text = value?.toString() ?: ''
        return "'${text.replace("'", "'\"'\"'")}'".toString()
    }
}
